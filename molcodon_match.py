#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from molcodon_trace import TraceToken, encode_with_trace

PHARMACOPHORE_MAP = {
    'OXN': 'HBA',
    'OXO': 'HBD',
    'CXN': 'POS',
    'CXS': 'POS',
    'CXO': 'NEG',
    'CXX': 'NEG',
}

_FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@dataclass
class TraceBundle:
    smiles: str
    tokens: List[str]
    trace: List[TraceToken]
    mol: object
    components: Dict[str, object]


@dataclass
class MatchResult:
    reference: TraceBundle
    hit: TraceBundle
    scores: Dict[str, object]
    matched_ref_atoms: Dict[int, str]
    matched_hit_atoms: Dict[int, str]
    matched_ref_bonds: Dict[int, str]
    matched_hit_bonds: Dict[int, str]
    matched_ref_token_spans: List[Tuple[int, int]]
    matched_hit_token_spans: List[Tuple[int, int]]
    matched_components: Dict[str, object]
    tanimoto: float


def calc_f1(match_count: int, query_total: int, target_total: int) -> float:
    if query_total == 0 and target_total == 0:
        return 100.0
    if query_total == 0 or target_total == 0 or match_count == 0:
        return 0.0
    precision = match_count / target_total
    recall = match_count / query_total
    return round((2.0 * precision * recall / (precision + recall)) * 100.0, 2)


def multiset_overlap(c1: Counter, c2: Counter) -> int:
    return sum((c1 & c2).values())


def sequence_similarity(seq1: Tuple[str, ...], seq2: Tuple[str, ...]) -> float:
    if not seq1 and not seq2:
        return 100.0
    if not seq1 or not seq2:
        return 0.0
    n = 3
    def grams(seq: Tuple[str, ...]) -> Counter:
        width = max(1, len(seq) - n + 1)
        return Counter(tuple(seq[i:i+n]) for i in range(width))
    g1 = grams(seq1)
    g2 = grams(seq2)
    inter = sum((g1 & g2).values())
    union = sum((g1 | g2).values())
    if union == 0:
        return 0.0
    return round((inter / union) * 100.0, 2)


def _token_span(indices: Iterable[int]) -> Optional[Tuple[int, int]]:
    idxs = sorted(set(indices))
    if not idxs:
        return None
    return (idxs[0], idxs[-1])


def _build_components(tokens: List[str], trace: List[TraceToken]) -> Dict[str, object]:
    rings: List[Tuple[str, ...]] = []
    branches: List[Tuple[str, ...]] = []
    ring_entries: List[Dict[str, object]] = []
    branch_entries: List[Dict[str, object]] = []
    ring_open_stack: List[int] = []
    branch_open_stack: List[int] = []

    groups: Dict[Tuple[str, int], List[TraceToken]] = defaultdict(list)
    atom_to_component: Dict[int, Set[Tuple[str, int]]] = defaultdict(set)
    bond_to_component: Dict[int, Set[Tuple[str, int]]] = defaultdict(set)
    pharmacophores: Counter = Counter()
    pharmacophore_entries: List[Dict[str, object]] = []
    bond_counter: Counter = Counter()
    bond_entries: List[Dict[str, object]] = []
    attachment_entries: List[Dict[str, object]] = []
    backbone_tokens: List[str] = []
    backbone_token_indices: List[int] = []
    backbone_atoms: Set[int] = set()
    backbone_bonds: Set[int] = set()

    for t in trace:
        if t.comp_type in ('ring', 'branch'):
            groups[(t.comp_type, t.comp_id)].append(t)
        if t.atom_idx is not None and t.comp_type in ('ring', 'branch'):
            atom_to_component[t.atom_idx].add((t.comp_type, t.comp_id))
        if t.bond_idx is not None and t.comp_type in ('ring', 'branch'):
            bond_to_component[t.bond_idx].add((t.comp_type, t.comp_id))

        if t.role == 'ring_open':
            ring_open_stack.append(t.idx)
        elif t.role == 'ring_close' and ring_open_stack:
            start = ring_open_stack.pop()
            ring_entries.append({'span': (start, t.idx), 'ring_id': t.comp_id})
        elif 'branch_open' in t.role:
            branch_open_stack.append(t.idx)
        elif 'branch_close' in t.role and branch_open_stack:
            start = branch_open_stack.pop()
            branch_entries.append({'span': (start, t.idx), 'branch_id': t.comp_id})

        if t.codon in PHARMACOPHORE_MAP and t.atom_idx is not None:
            label = PHARMACOPHORE_MAP[t.codon]
            context = 'backbone'
            if t.comp_type == 'ring':
                context = 'ring'
            elif t.comp_type == 'branch':
                context = 'branch'
            # Pharmacophore signature includes context type only.
            # The specific component identity (which ring / which branch)
            # is validated through the Ring F1 and Branch F1 scores;
            # pharmacophore F1 measures whether the same feature types
            # exist in the same structural context categories.
            sig = (label, context)
            pharmacophores[sig] += 1
            pharmacophore_entries.append({
                'signature': sig,
                'atom_idx': t.atom_idx,
                'token_idx': t.idx,
                'context': context,
            })

        if t.role == 'bond' and t.bond_idx is not None:
            bond_counter[t.codon] += 1
            bond_entries.append({
                'codon': t.codon,
                'bond_idx': t.bond_idx,
                'token_idx': t.idx,
                'context': t.comp_type,
                'comp_id': t.comp_id,
            })

        in_component = False
        if t.atom_idx is not None and atom_to_component.get(t.atom_idx):
            in_component = True
        if t.bond_idx is not None and bond_to_component.get(t.bond_idx):
            in_component = True
        if not in_component and t.codon not in ('SCC', 'SSS') and t.role not in {'start', 'end'}:
            backbone_tokens.append(t.codon)
            backbone_token_indices.append(t.idx)
            if t.atom_idx is not None:
                backbone_atoms.add(t.atom_idx)
            if t.bond_idx is not None:
                backbone_bonds.add(t.bond_idx)

    ATOM_CODONS = {'CCC', 'CCN', 'CCO', 'CCS', 'CNC', 'CNN', 'CNO', 'CNS', 'COC', 'CON'}
    BOND_CODONS = {'NCC', 'NCN', 'NCO', 'NCS'}

    def _normalize_ring_signature(sig: Tuple[str, ...]) -> Tuple[str, ...]:
        """Normalize ring signatures for topology-only comparison.

        Keeps ONLY atom codons and bond codons. Strips everything else:
        - Charge annotations (CCX, CXN, ...) → Pharmacophore F1 handles charged atoms
        - Pharmacophore annotations (OXN, OXO) → Pharmacophore F1 handles HBA/HBD
        - Bond annotations (NXO, NXS, ...) → aromatic info preserved via bond codon normalization
        - Fusion markers (OCC, SOC, CSO, ...) → fused atoms compared by atom type only
        - Stereo annotations (SXN, SXO) → separate concern

        Before stripping, aromatic bonds (Kekulized bond + NXO) are normalized to NCS.
        """
        # Step 1: normalize aromatic bonds before stripping annotations
        temp = list(sig)
        i = 0
        while i < len(temp) - 1:
            if temp[i] in BOND_CODONS and temp[i + 1] == 'NXO':
                temp[i] = 'NCS'
            i += 1
        # Step 2: keep only atom and bond codons (pure topology)
        topo = tuple(c for c in temp if c in ATOM_CODONS or c in BOND_CODONS)
        # Step 3: canonical rotation — the same ring entered from different
        # atoms produces different linear sequences.  We pick the
        # lexicographically smallest rotation so that the signature is
        # entry-point-invariant.  Rotation unit = 2 codons (atom + bond).
        n = len(topo)
        if n < 4:
            return topo
        step = 2
        best = topo
        for start in range(step, n, step):
            rotated = topo[start:] + topo[:start]
            if rotated < best:
                best = rotated
        return best
        return tuple(result)

    def build_component_entries(kind: str, items: Dict[Tuple[str, int], List[TraceToken]]) -> List[Dict[str, object]]:
        out: List[Dict[str, object]] = []
        for (ctype, cid), toks in items.items():
            if ctype != kind:
                continue
            payload = tuple(tok.codon for tok in toks if tok.role not in {'ring_open', 'ring_close', 'branch_open', 'branch_close', 'branch_replay_open', 'branch_replay_close', 'branch_empty_open', 'branch_empty_close'})
            if kind == 'ring':
                payload = _normalize_ring_signature(payload)
            atoms = sorted({tok.atom_idx for tok in toks if tok.atom_idx is not None})
            bonds = sorted({tok.bond_idx for tok in toks if tok.bond_idx is not None})
            span = _token_span(tok.idx for tok in toks)
            out.append({
                'id': cid,
                'signature': payload,
                'atoms': atoms,
                'bonds': bonds,
                'span': span,
                'token_indices': [tok.idx for tok in toks],
            })
        out.sort(key=lambda x: x['id'])
        return out

    ring_components = build_component_entries('ring', groups)
    branch_components = build_component_entries('branch', groups)

    component_by_atom: Dict[int, Tuple[str, int]] = {}
    for entry in ring_components:
        for atom_idx in entry['atoms']:
            component_by_atom.setdefault(atom_idx, ('ring', entry['id']))
    for entry in branch_components:
        for atom_idx in entry['atoms']:
            component_by_atom.setdefault(atom_idx, ('branch', entry['id']))

    entry_lookup = {
        ('ring', entry['id']): entry for entry in ring_components
    }
    entry_lookup.update({('branch', entry['id']): entry for entry in branch_components})

    # ── Phase 1: find each branch's parent via trace ────────────
    branch_by_id: Dict[int, Dict] = {b['id']: b for b in branch_components}
    branch_parent: Dict[int, Dict] = {}  # branch_id → {parent_type, parent_id, attach_atom, ...}

    for branch in branch_components:
        first_bond_codon = None
        for tok in trace:
            if tok.comp_type == 'branch' and tok.comp_id == branch['id']:
                if first_bond_codon is None and tok.role == 'bond':
                    first_bond_codon = tok.codon

        parent_type = 'backbone'
        parent_id = None
        attach_atom = None
        parent_signature = ('backbone',)
        branch_tokens = groups.get(('branch', branch['id']), [])
        if branch_tokens:
            first_idx = branch_tokens[0].idx
            if first_idx > 0:
                opener = trace[first_idx - 1]
                if opener.atom_idx is not None:
                    attach_atom = opener.atom_idx
                    if attach_atom in component_by_atom:
                        parent_type, parent_id = component_by_atom[attach_atom]
                        parent_entry = entry_lookup.get((parent_type, parent_id))
                        if parent_entry:
                            parent_signature = parent_entry['signature']

        branch['parent_type'] = parent_type
        branch_parent[branch['id']] = {
            'parent_type': parent_type,
            'parent_id': parent_id,
            'attach_atom': attach_atom,
            'parent_signature': parent_signature,
            'first_bond_codon': first_bond_codon,
        }

    # ── Phase 2: build parent→children tree ───────────────────
    branch_children: Dict[int, List[int]] = defaultdict(list)
    for bid, pinfo in branch_parent.items():
        if pinfo['parent_type'] == 'branch' and pinfo['parent_id'] is not None:
            branch_children[pinfo['parent_id']].append(bid)

    # ── Phase 2b: detect destination ring for each branch ─────
    # If a branch leads to a ring at its far end, record that
    # ring's normalized signature as the "destination".
    ring_entry_by_id = {r['id']: r for r in ring_components}
    for branch in branch_components:
        branch_toks = groups.get(('branch', branch['id']), [])
        branch['destination'] = None
        if branch_toks:
            last_idx = branch_toks[-1].idx
            # Scan forward from the branch's last token to find a ring entry
            for offset in range(1, 4):
                nidx = last_idx + offset
                if nidx >= len(trace):
                    break
                nt = trace[nidx]
                if nt.comp_type == 'ring' and nt.comp_id is not None:
                    re = ring_entry_by_id.get(nt.comp_id)
                    if re:
                        branch['destination'] = re['signature']
                    break

    # ── Phase 3: recursive fragment signature ─────────────────
    # Fragment = (own_signature, child_fragments, destination_ring_sig_or_None)
    def _fragment_sig(bid: int) -> tuple:
        branch = branch_by_id.get(bid)
        if not branch:
            return ()
        children = branch_children.get(bid, [])
        children.sort(key=lambda c: branch_by_id.get(c, {}).get('signature', ()))
        child_sigs = tuple(_fragment_sig(c) for c in children)
        return (branch['signature'], child_sigs, branch.get('destination'))

    # ── Phase 4: set match_key and build attachment_entries ───
    for branch in branch_components:
        pinfo = branch_parent[branch['id']]
        frag = _fragment_sig(branch['id'])

        # Branch match_key: parent_type + OWN content only (children matched separately)
        branch['match_key'] = (pinfo['parent_type'],) + branch['signature']

        attachment_entries.append({
            'branch_id': branch['id'],
            'parent_type': pinfo['parent_type'],
            'parent_id': pinfo['parent_id'],
            'attach_atom': pinfo['attach_atom'],
            'parent_signature': pinfo['parent_signature'],
            'fragment': frag,
            'signature': (pinfo['parent_type'], pinfo['parent_signature'], frag),
            'branch_signature': branch['signature'],
        })

    rings = Counter(entry['signature'] for entry in ring_components)
    branches = Counter(entry.get('match_key', entry['signature']) for entry in branch_components)
    attachments = Counter(entry['signature'] for entry in attachment_entries)

    return {
        'rings': rings,
        'ring_entries': ring_components,
        'branches': branches,
        'branch_entries': branch_components,
        'bond_types': bond_counter,
        'bond_entries': bond_entries,
        'attachments': attachments,
        'attachment_entries': attachment_entries,
        'backbone': tuple(backbone_tokens),
        'backbone_entries': [{
            'signature': tuple(backbone_tokens),
            'span': _token_span(backbone_token_indices),
            'token_indices': backbone_token_indices,
            'atoms': sorted(backbone_atoms),
            'bonds': sorted(backbone_bonds),
        }],
        'pharmacophores': pharmacophores,
        'pharmacophore_entries': pharmacophore_entries,
        'ring_spans': [e['span'] for e in ring_entries if e.get('span')],
        'branch_spans': [e['span'] for e in branch_entries if e.get('span')],
    }


def build_trace_bundle(smiles: str) -> TraceBundle:
    tokens, trace, mol = encode_with_trace(smiles)
    return TraceBundle(
        smiles=smiles,
        tokens=tokens,
        trace=trace,
        mol=mol,
        components=_build_components(tokens, trace),
    )


def _match_key(entry: Dict[str, object]) -> tuple:
    return entry.get('match_key', entry['signature'])


def _match_component_entries(ref_entries: List[Dict[str, object]], hit_entries: List[Dict[str, object]]) -> Tuple[List[Tuple[Dict[str, object], Dict[str, object]]], int]:
    by_sig: Dict[tuple, List[Dict[str, object]]] = defaultdict(list)
    for entry in hit_entries:
        by_sig[_match_key(entry)].append(entry)
    pairs: List[Tuple[Dict[str, object], Dict[str, object]]] = []
    match_count = 0
    for ref_entry in ref_entries:
        bucket = by_sig.get(_match_key(ref_entry))
        if bucket:
            hit_entry = bucket.pop(0)
            pairs.append((ref_entry, hit_entry))
            match_count += 1
    return pairs, match_count


def _match_rings_contextual(
    ref_entries: List[Dict], hit_entries: List[Dict],
    ref_attachments: List[Dict], hit_attachments: List[Dict],
) -> Tuple[List[Tuple[Dict, Dict]], int]:
    """Match rings by topology, using branch decoration as tie-breaker.

    When multiple hit rings have the same topology as a ref ring,
    pick the one whose attached branches are most similar.
    """
    # Build ring_id → set of branch signatures attached to it
    def _ring_branches(att_list: List[Dict]) -> Dict[int, Set[tuple]]:
        out: Dict[int, Set[tuple]] = defaultdict(set)
        for a in att_list:
            if a['parent_type'] == 'ring' and a['parent_id'] is not None:
                out[a['parent_id']].add(a.get('branch_signature', ()))
        return out

    ref_rb = _ring_branches(ref_attachments)
    hit_rb = _ring_branches(hit_attachments)

    by_sig: Dict[tuple, List[Dict]] = defaultdict(list)
    for entry in hit_entries:
        by_sig[entry['signature']].append(entry)

    pairs: List[Tuple[Dict, Dict]] = []
    match_count = 0
    for ref_entry in ref_entries:
        bucket = by_sig.get(ref_entry['signature'])
        if not bucket:
            continue

        if len(bucket) == 1:
            pairs.append((ref_entry, bucket.pop(0)))
            match_count += 1
            continue

        # Multiple candidates → pick by branch decoration overlap
        r_branches = ref_rb.get(ref_entry['id'], set())
        best_score = -1.0
        best_idx = 0
        for i, hit_entry in enumerate(bucket):
            h_branches = hit_rb.get(hit_entry['id'], set())
            if not r_branches and not h_branches:
                score = 1.0
            elif not r_branches or not h_branches:
                score = 0.0
            else:
                overlap = len(r_branches & h_branches)
                total = len(r_branches | h_branches)
                score = overlap / total
            if score > best_score:
                best_score = score
                best_idx = i

        pairs.append((ref_entry, bucket.pop(best_idx)))
        match_count += 1

    return pairs, match_count


def _match_branches_contextual(
    ref_entries: List[Dict], hit_entries: List[Dict],
    ref_attachments: List[Dict], hit_attachments: List[Dict],
    ring_pairs: List[Tuple[Dict, Dict]],
) -> Tuple[List[Tuple[Dict, Dict]], int]:
    """Match branches by content, using ring-pair context and attachment as tie-breakers.

    Prefers matching branches that sit on the same matched ring pair.
    """
    ref_att = {a['branch_id']: a for a in ref_attachments}
    hit_att = {a['branch_id']: a for a in hit_attachments}

    # Build ring ID mapping from ring_pairs: ref_ring_id → hit_ring_id
    ring_map: Dict[int, int] = {}
    for r_ring, h_ring in ring_pairs:
        ring_map[r_ring['id']] = h_ring['id']

    by_key: Dict[tuple, List[Dict]] = defaultdict(list)
    for entry in hit_entries:
        by_key[_match_key(entry)].append(entry)

    pairs: List[Tuple[Dict, Dict]] = []
    match_count = 0
    for ref_entry in ref_entries:
        bucket = by_key.get(_match_key(ref_entry))
        if not bucket:
            continue

        if len(bucket) == 1:
            pairs.append((ref_entry, bucket.pop(0)))
            match_count += 1
            continue

        # Multiple candidates → score each
        ra = ref_att.get(ref_entry['id'])
        r_frag = ra.get('fragment', ()) if ra else ()
        r_psig = ra.get('parent_signature', ()) if ra else ()
        r_parent_id = ra.get('parent_id') if ra else None

        best_score = -1.0
        best_idx = 0
        for i, hit_entry in enumerate(bucket):
            ha = hit_att.get(hit_entry['id'])
            h_frag = ha.get('fragment', ()) if ha else ()
            h_psig = ha.get('parent_signature', ()) if ha else ()
            h_parent_id = ha.get('parent_id') if ha else None

            score = _fragment_similarity(r_frag, h_frag)
            if r_psig == h_psig:
                score += 1.0
            # Strong bonus: both branches sit on matched ring pair
            if r_parent_id is not None and ring_map.get(r_parent_id) == h_parent_id:
                score += 2.0
            if score > best_score:
                best_score = score
                best_idx = i

        pairs.append((ref_entry, bucket.pop(best_idx)))
        match_count += 1

    return pairs, match_count


def _pair_pharmacophores(ref_entries: List[Dict[str, object]], hit_entries: List[Dict[str, object]]) -> List[Tuple[Dict[str, object], Dict[str, object]]]:
    by_sig: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for entry in hit_entries:
        by_sig[entry['signature']].append(entry)
    pairs = []
    for ref_entry in ref_entries:
        bucket = by_sig.get(ref_entry['signature'])
        if bucket:
            pairs.append((ref_entry, bucket.pop(0)))
    return pairs


def _fragment_similarity(frag_a: tuple, frag_b: tuple) -> float:
    """Recursive similarity between two fragment trees (0.0 – 1.0).

    A fragment is ``(own_signature, (child_fragments...), destination_or_None)``.
    - Own content must match for any credit (base 50%).
    - Children and destination contribute the remaining 50%.
    """
    if not frag_a and not frag_b:
        return 1.0
    if not frag_a or not frag_b:
        return 0.0

    sig_a = frag_a[0]
    sig_b = frag_b[0]
    children_a = frag_a[1] if len(frag_a) > 1 else ()
    children_b = frag_b[1] if len(frag_b) > 1 else ()
    dest_a = frag_a[2] if len(frag_a) > 2 else None
    dest_b = frag_b[2] if len(frag_b) > 2 else None

    if sig_a != sig_b:
        return 0.0

    # Own content matches → collect sub-scores for children + destination
    sub_scores: List[float] = []

    # Children comparison
    if children_a or children_b:
        if not children_a or not children_b:
            sub_scores.append(0.0)
        else:
            remaining = list(children_b)
            child_scores: List[float] = []
            for ca in children_a:
                best_score = 0.0
                best_idx = -1
                for i, cb in enumerate(remaining):
                    s = _fragment_similarity(ca, cb)
                    if s > best_score:
                        best_score = s
                        best_idx = i
                if best_idx >= 0 and best_score > 0:
                    remaining.pop(best_idx)
                child_scores.append(best_score)
            total = max(len(children_a), len(children_b))
            sub_scores.append(sum(child_scores) / total)

    # Destination ring comparison
    if dest_a is not None or dest_b is not None:
        sub_scores.append(1.0 if dest_a == dest_b else 0.0)

    if not sub_scores:
        return 1.0  # leaf node, no children, no destination → perfect

    return 0.5 + 0.5 * (sum(sub_scores) / len(sub_scores))


def _match_attachments_graded(ref_entries: List[Dict], hit_entries: List[Dict]) -> Tuple[float, List[Tuple[Dict, Dict, float]]]:
    """Match attachments with graded fragment similarity.

    Returns (f1_score, pairs) where each pair is (ref_entry, hit_entry, similarity).
    """
    if not ref_entries and not hit_entries:
        return 100.0, []
    if not ref_entries or not hit_entries:
        return 0.0, []

    ref_total = len(ref_entries)
    hit_total = len(hit_entries)

    hit_available = list(range(len(hit_entries)))
    total_sim = 0.0
    pairs: List[Tuple[Dict, Dict, float]] = []

    for ref_entry in ref_entries:
        best_sim = 0.0
        best_idx = -1

        for i in hit_available:
            hit_entry = hit_entries[i]
            if ref_entry['parent_type'] != hit_entry['parent_type']:
                continue

            frag_sim = _fragment_similarity(
                ref_entry.get('fragment', ()),
                hit_entry.get('fragment', ()),
            )
            if frag_sim == 0.0:
                continue

            if ref_entry.get('parent_signature') == hit_entry.get('parent_signature'):
                sim = frag_sim
            else:
                sim = frag_sim * 0.5

            if sim > best_sim:
                best_sim = sim
                best_idx = i

        if best_idx >= 0:
            hit_available.remove(best_idx)
            pairs.append((ref_entry, hit_entries[best_idx], best_sim))
        total_sim += best_sim

    precision = total_sim / hit_total
    recall = total_sim / ref_total
    if precision + recall == 0:
        return 0.0, pairs
    return round(200.0 * precision * recall / (precision + recall), 2), pairs


def tanimoto_score(smiles_a: str, smiles_b: str) -> float:
    try:
        ma = Chem.MolFromSmiles(smiles_a)
        mb = Chem.MolFromSmiles(smiles_b)
        if ma is None or mb is None:
            return 0.0
        fp_a = _FPGEN.GetFingerprint(ma)
        fp_b = _FPGEN.GetFingerprint(mb)
        return round(100.0 * DataStructs.TanimotoSimilarity(fp_a, fp_b), 2)
    except Exception:
        return 0.0


def match_smiles(reference_smiles: str, hit_smiles: str) -> MatchResult:
    reference = build_trace_bundle(reference_smiles)
    hit = build_trace_bundle(hit_smiles)
    ref_comp = reference.components
    hit_comp = hit.components

    ring_pairs, ring_match = _match_rings_contextual(
        ref_comp['ring_entries'], hit_comp['ring_entries'],
        ref_comp['attachment_entries'], hit_comp['attachment_entries'],
    )
    branch_pairs, branch_match = _match_branches_contextual(
        ref_comp['branch_entries'], hit_comp['branch_entries'],
        ref_comp['attachment_entries'], hit_comp['attachment_entries'],
        ring_pairs,
    )
    pharm_pairs = _pair_pharmacophores(ref_comp['pharmacophore_entries'], hit_comp['pharmacophore_entries'])

    bond_match = multiset_overlap(ref_comp['bond_types'], hit_comp['bond_types'])
    bond_q_tot = sum(ref_comp['bond_types'].values())
    bond_t_tot = sum(hit_comp['bond_types'].values())

    pharm_match = multiset_overlap(ref_comp['pharmacophores'], hit_comp['pharmacophores'])
    pharm_q_tot = sum(ref_comp['pharmacophores'].values())
    pharm_t_tot = sum(hit_comp['pharmacophores'].values())

    ring_q_tot = sum(ref_comp['rings'].values())
    ring_t_tot = sum(hit_comp['rings'].values())
    branch_q_tot = sum(ref_comp['branches'].values())
    branch_t_tot = sum(hit_comp['branches'].values())

    attachment_q_tot = len(ref_comp['attachment_entries'])
    attachment_t_tot = len(hit_comp['attachment_entries'])

    ring_f1 = calc_f1(ring_match, ring_q_tot, ring_t_tot)
    branch_f1 = calc_f1(branch_match, branch_q_tot, branch_t_tot)
    attachment_f1, attachment_graded_pairs = _match_attachments_graded(ref_comp['attachment_entries'], hit_comp['attachment_entries'])
    bond_f1 = calc_f1(bond_match, bond_q_tot, bond_t_tot)
    pharm_f1 = calc_f1(pharm_match, pharm_q_tot, pharm_t_tot)
    backbone_score = sequence_similarity(ref_comp['backbone'], hit_comp['backbone'])

    overall = round((ring_f1 + branch_f1 + 2.0 * attachment_f1 + pharm_f1) / 5.0, 2)

    ref_atoms: Dict[int, str] = {}
    hit_atoms: Dict[int, str] = {}
    ref_bonds: Dict[int, str] = {}
    hit_bonds: Dict[int, str] = {}
    ref_spans: List[Tuple[int, int]] = []
    hit_spans: List[Tuple[int, int]] = []

    def add_match(atoms_d, bonds_d, atoms, bonds, role):
        for a in atoms:
            if a not in atoms_d: atoms_d[a] = role
        for b in bonds:
            if b not in bonds_d: bonds_d[b] = role

    for ref_entry, hit_entry in ring_pairs:
        add_match(ref_atoms, ref_bonds, ref_entry['atoms'], ref_entry['bonds'], 'ring')
        add_match(hit_atoms, hit_bonds, hit_entry['atoms'], hit_entry['bonds'], 'ring')
        if ref_entry.get('span'): ref_spans.append(ref_entry['span'])
        if hit_entry.get('span'): hit_spans.append(hit_entry['span'])

    for ref_entry, hit_entry in branch_pairs:
        add_match(ref_atoms, ref_bonds, ref_entry['atoms'], ref_entry['bonds'], 'branch')
        add_match(hit_atoms, hit_bonds, hit_entry['atoms'], hit_entry['bonds'], 'branch')
        if ref_entry.get('span'): ref_spans.append(ref_entry['span'])
        if hit_entry.get('span'): hit_spans.append(hit_entry['span'])

    for ref_entry, hit_entry in pharm_pairs:
        add_match(ref_atoms, ref_bonds, [ref_entry['atom_idx']], [], 'pharm')
        add_match(hit_atoms, hit_bonds, [hit_entry['atom_idx']], [], 'pharm')
        ref_spans.append((ref_entry['token_idx'], ref_entry['token_idx']))
        hit_spans.append((hit_entry['token_idx'], hit_entry['token_idx']))

    ref_backbone = ref_comp['backbone_entries'][0]
    hit_backbone = hit_comp['backbone_entries'][0]

    ref_bond_entries_by_codon: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    hit_bond_entries_by_codon: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for entry in ref_comp['bond_entries']:
        ref_bond_entries_by_codon[entry['codon']].append(entry)
    for entry in hit_comp['bond_entries']:
        hit_bond_entries_by_codon[entry['codon']].append(entry)
    for codon, ref_entries in ref_bond_entries_by_codon.items():
        hit_entries = hit_bond_entries_by_codon.get(codon, [])
        for ref_entry, hit_entry in zip(ref_entries, hit_entries):
            add_match(ref_atoms, ref_bonds, [], [ref_entry['bond_idx']], 'bond')
            add_match(hit_atoms, hit_bonds, [], [hit_entry['bond_idx']], 'bond')
            ref_spans.append((ref_entry['token_idx'], ref_entry['token_idx']))
            hit_spans.append((hit_entry['token_idx'], hit_entry['token_idx']))

    scores = {
        'overall': overall,
        'ring_f1': ring_f1,
        'ring_match': f'{ring_match}/{ring_q_tot}',
        'branch_f1': branch_f1,
        'branch_match': f'{branch_match}/{branch_q_tot}',
        'attachment_f1': attachment_f1,
        'attachment_match': f'{len(attachment_graded_pairs)}/{max(attachment_q_tot, attachment_t_tot)}',
        'bond_f1': bond_f1,
        'bond_match': f'{bond_match}/{bond_q_tot}',
        'pharmacophore_f1': pharm_f1,
        'pharmacophore_match': f'{pharm_match}/{pharm_q_tot}',
        'backbone': backbone_score,
    }

    matched_components = {
        'rings': [{'ref_id': a['id'], 'hit_id': b['id'], 'ref_atoms': a['atoms'], 'hit_atoms': b['atoms'], 'ref_bonds': a['bonds'], 'hit_bonds': b['bonds'], 'signature': list(a['signature'])} for a, b in ring_pairs],
        'branches': [{'ref_id': a['id'], 'hit_id': b['id'], 'ref_atoms': a['atoms'], 'hit_atoms': b['atoms'], 'ref_bonds': a['bonds'], 'hit_bonds': b['bonds'], 'parent_type': a.get('parent_type', ''), 'signature': list(a['signature'])} for a, b in branch_pairs],
        'attachments': [{'ref_branch_id': a['branch_id'], 'hit_branch_id': b['branch_id'], 'parent_type': a['parent_type'], 'ref_attach_atom': a.get('attach_atom'), 'hit_attach_atom': b.get('attach_atom'), 'similarity': round(sim, 2)} for a, b, sim in attachment_graded_pairs],
        'pharmacophores': [{'label': a['signature'][0], 'context': a['signature'][1], 'ref_atom': a['atom_idx'], 'hit_atom': b['atom_idx']} for a, b in pharm_pairs],
    }

    return MatchResult(
        reference=reference,
        hit=hit,
        scores=scores,
        matched_ref_atoms=ref_atoms,
        matched_hit_atoms=hit_atoms,
        matched_ref_bonds=ref_bonds,
        matched_hit_bonds=hit_bonds,
        matched_ref_token_spans=sorted(set(ref_spans)),
        matched_hit_token_spans=sorted(set(hit_spans)),
        matched_components=matched_components,
        tanimoto=tanimoto_score(reference_smiles, hit_smiles),
    )
