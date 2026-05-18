import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Optional, Set

from rdkit import Chem
from molcodon_encoder import (
    ATOM_TO_CODON, BOND_TYPE_TO_CODON, BRANCH_OPEN, BRANCH_CLOSE,
    RING_OPEN, RING_CLOSE, START_CODON, END_CODON, FUSION_MARKER,
    POSITION_CODONS, RING_REF_CODONS, EncodeError,
    get_atom_annotations, get_bond_annotations,
    atom_order_key, neighbor_order_key
)

@dataclass
class TraceToken:
    idx: int
    codon: str
    role: str
    atom_idx: Optional[int] = None
    bond_idx: Optional[int] = None
    comp_type: str = 'backbone'
    comp_id: Optional[int] = None

class Tracer:
    def __init__(self):
        self.tokens: List[str] = []
        self.trace: List[TraceToken] = []
        self.comp_type = 'backbone'
        self.comp_id = 0
        self.comp_stack = []
        
    def push_comp(self, ctype, cid):
        self.comp_stack.append((self.comp_type, self.comp_id))
        self.comp_type = ctype
        self.comp_id = cid
        
    def pop_comp(self):
        self.comp_type, self.comp_id = self.comp_stack.pop()

    def emit(self, codon, role, atom_idx=None, bond_idx=None):
        idx = len(self.tokens)
        self.tokens.append(codon)
        self.trace.append(TraceToken(
            idx=idx, codon=codon, role=role,
            atom_idx=atom_idx, bond_idx=bond_idx,
            comp_type=self.comp_type, comp_id=self.comp_id
        ))

def encode_with_trace(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: raise EncodeError(f"Invalid SMILES: {smiles}")
    if mol.GetNumAtoms() == 0: raise EncodeError("Empty molecule")
    atom_h_counts = {a.GetIdx(): a.GetTotalNumHs() for a in mol.GetAtoms()}
    aromatic_bonds = {b.GetIdx() for b in mol.GetBonds() if b.GetIsAromatic()}
    Chem.Kekulize(mol, clearAromaticFlags=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    canon = list(Chem.CanonicalRankAtoms(mol))
    all_rings = mol.GetRingInfo().AtomRings()
    a2r = defaultdict(list)
    for ri, ring in enumerate(all_rings):
        for ai in ring: a2r[ai].append(ri)

    adj = defaultdict(list)
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj[a1].append((a2, b)); adj[a2].append((a1, b))

    tracer = Tracer()
    tracer.emit(START_CODON, 'start')

    vis_a, vis_r = set(), set()
    _slot_next = [0]; _slot_free = []
    def alloc():
        if _slot_free: return _slot_free.pop(0)
        s = _slot_next[0]
        if s >= len(RING_OPEN): raise EncodeError("Ring label capacity exceeded")
        _slot_next[0] += 1
        return s
    def free(s):
        if s not in _slot_free: _slot_free.append(s)

    completed_rings = []
    def uring(ai):
        for ri in a2r.get(ai, []):
            if ri not in vis_r: return ri
        return None

    def neighbor_key(cur, nb, bond):
        return neighbor_order_key(mol, canon, nb, bond)
    start = min(range(mol.GetNumAtoms()), key=lambda ai: atom_order_key(mol, canon, ai))

    def rorder(ri, entry):
        rs = set(all_rings[ri])
        if entry not in rs: return list(rs)
        ring_adj = {ai: [(n, b) for n, b in adj[ai] if n in rs] for ai in rs}
        next_candidates = ring_adj[entry]
        if not next_candidates: return [entry]
        next_candidates.sort(key=lambda x: neighbor_order_key(mol, canon, x[0], x[1]))
        first_next = next_candidates[0][0]
        o, cur, prev = [entry, first_next], first_next, entry
        while len(o) < len(rs):
            nb = [(n, b) for n, b in ring_adj[cur] if n != prev and n not in o]
            if not nb: break
            nb.sort(key=lambda x: neighbor_order_key(mol, canon, x[0], x[1]))
            o.append(nb[0][0])
            prev, cur = cur, nb[0][0]
        return o

    def emit_bond(bond):
        if not bond: return
        bidx = bond.GetIdx()
        bt = bond.GetBondType()
        tracer.emit(BOND_TYPE_TO_CODON[bt], 'bond', bond_idx=bidx)
        for ann in get_bond_annotations(bond, aromatic_bonds, mol=mol, canon=canon):
            tracer.emit(ann, 'bond_annotation', bond_idx=bidx)

    def emit_atom(ai):
        atom = mol.GetAtomWithIdx(ai)
        sym = atom.GetSymbol()
        tracer.emit(ATOM_TO_CODON[sym], 'atom', atom_idx=ai)
        for ann in get_atom_annotations(atom, atom_h_counts):
            tracer.emit(ann, 'atom_annotation', atom_idx=ai)

    branch_counter = [0]
    
    def traverse(ai, from_bond=None):
        if ai in vis_a: return
        ri = uring(ai)
        if ri is not None:
            enc_ring(ri, ai, from_bond)
            return
        vis_a.add(ai)
        if from_bond is not None: emit_bond(from_bond)
        emit_atom(ai)

        unvis = [(n, b) for n, b in adj[ai] if n not in vis_a]
        if len(unvis) == 1:
            traverse(unvis[0][0], unvis[0][1])
        elif len(unvis) > 1:
            unvis.sort(key=lambda x: neighbor_key(ai, x[0], x[1]))
            for i, (n, b) in enumerate(unvis):
                if i < len(unvis) - 1:
                    d = i % len(BRANCH_OPEN)
                    b_id = branch_counter[0]; branch_counter[0] += 1
                    tracer.emit(BRANCH_OPEN[d], 'branch_open', atom_idx=ai)
                    tracer.push_comp('branch', b_id)
                    traverse(n, b)
                    tracer.pop_comp()
                    tracer.emit(BRANCH_CLOSE[d], 'branch_close', atom_idx=ai)
                else:
                    traverse(n, b)

    def enc_ring(ri, entry, from_bond):
        if ri in vis_r: return
        vis_r.add(ri)
        if from_bond is not None: emit_bond(from_bond)
        tracer.push_comp('ring', ri)
        slot = alloc()
        tracer.emit(RING_OPEN[slot], 'ring_open', atom_idx=entry)

        ring_atoms = rorder(ri, entry)
        ring_set = set(ring_atoms)
        ring_atom_list = []
        deferred = []
        seen_fused_targets = set()
        
        prev = None
        for i, a in enumerate(ring_atoms):
            already_visited = a in vis_a
            vis_a.add(a)
            if i > 0 and prev is not None:
                b = mol.GetBondBetweenAtoms(prev, a)
                if b: emit_bond(b)
            
            if already_visited:
                tracer.emit(FUSION_MARKER, 'fusion', atom_idx=a)
                fusion_written = False
                for cri, cring in enumerate(completed_rings):
                    if a in cring:
                        pos = cring.index(a)
                        tracer.emit(RING_REF_CODONS[cri], 'fusion_ref', atom_idx=a)
                        tracer.emit(POSITION_CODONS[pos], 'fusion_pos', atom_idx=a)
                        fusion_written = True
                        break
                atom = mol.GetAtomWithIdx(a)
                tracer.emit(ATOM_TO_CODON[atom.GetSymbol()], 'atom', atom_idx=a)
            else:
                emit_atom(a)
            ring_atom_list.append(a)
            
            if not already_visited:
                for nb, nb_bond in adj[a]:
                    if nb in vis_a or nb in ring_set: continue
                    target_ring = uring(nb)
                    if target_ring is not None:
                        if target_ring not in seen_fused_targets:
                            seen_fused_targets.add(target_ring)
                            fused_label_slot = alloc()
                            deferred.append(('fused', fused_label_slot, a, target_ring, nb, nb_bond))
                    else:
                        d = branch_counter[0] % len(BRANCH_OPEN)
                        b_id = branch_counter[0]; branch_counter[0] += 1
                        tracer.emit(BRANCH_OPEN[d], 'branch_open', atom_idx=a)
                        deferred.append(('branch', d, a, nb, nb_bond, b_id))
            prev = a

        closing_bond = mol.GetBondBetweenAtoms(ring_atoms[-1], ring_atoms[0])
        if closing_bond: emit_bond(closing_bond)
        tracer.emit(RING_CLOSE[slot], 'ring_close')
        free(slot)
        tracer.pop_comp()
        completed_rings.append(ring_atom_list)
        
        for item in deferred:
            if item[0] == 'branch':
                _, br_idx, from_atom, nb, nb_bond, b_id = item
                if nb not in vis_a:
                    tracer.emit(BRANCH_OPEN[br_idx], 'branch_replay_open', atom_idx=from_atom)
                    tracer.push_comp('branch', b_id)
                    traverse(nb, nb_bond)
                    tracer.pop_comp()
                    tracer.emit(BRANCH_CLOSE[br_idx], 'branch_replay_close', atom_idx=from_atom)
                else:
                    tracer.emit(BRANCH_OPEN[br_idx], 'branch_empty_open', atom_idx=from_atom)
                    tracer.emit(BRANCH_CLOSE[br_idx], 'branch_empty_close', atom_idx=from_atom)
            elif item[0] == 'fused':
                _, fused_label_slot, from_atom, target_ring, nb, nb_bond = item
                if target_ring not in vis_r:
                    parent_ring_idx = None
                    pos = None
                    for cri, cring in enumerate(completed_rings):
                        if from_atom in cring:
                            parent_ring_idx = cri; pos = cring.index(from_atom); break
                    tracer.emit(RING_REF_CODONS[parent_ring_idx], 'fused_replay_ref', atom_idx=from_atom)
                    tracer.emit(POSITION_CODONS[pos], 'fused_replay_pos', atom_idx=from_atom)
                    _slot_free.insert(0, fused_label_slot)
                    actual_bond = mol.GetBondBetweenAtoms(from_atom, nb)
                    enc_ring(target_ring, nb, actual_bond)
                else:
                    free(fused_label_slot)

    traverse(start)
    tracer.emit(END_CODON, 'end')
    return tracer.tokens, tracer.trace, mol

if __name__ == '__main__':
    tokens, trace, m = encode_with_trace('c1ccccc1O')
    print(" ".join(tokens))
    for t in trace:
        print(f"{t.idx:02d} {t.codon} role={t.role:15} comp={t.comp_type}:{t.comp_id} a={t.atom_idx} b={t.bond_idx}")
