#!/usr/bin/env python3
import argparse
import sys
from rdkit import Chem
from collections import defaultdict


class EncodeError(Exception):
    pass

ATOM_TO_CODON = {
    'C': 'CCC', 'N': 'CCN', 'O': 'CCO', 'S': 'CCS',
    'F': 'CNC', 'Cl': 'CNN', 'Br': 'CNO', 'I': 'CNS',
    'P': 'COC', 'B': 'CON',
}
CODON_TO_ATOM = {v: k for k, v in ATOM_TO_CODON.items()}

BOND_TYPE_TO_CODON = {
    Chem.BondType.SINGLE: 'NCC', Chem.BondType.DOUBLE: 'NCN',
    Chem.BondType.TRIPLE: 'NCO', Chem.BondType.AROMATIC: 'NCS',
}

BRANCH_OPEN  = ['NNC', 'NOC', 'NOS', 'NSC']
BRANCH_CLOSE = ['NNN', 'NON', 'NOO', 'NSN']

RING_OPEN  = ['NNO', 'NSO', 'OSO', 'OCN', 'OCO', 'OCS', 'ONC', 'ONN']
RING_CLOSE = ['NNS', 'NSS', 'OSS', 'ONO', 'ONS', 'OSC', 'OSN', 'SCN']

START_CODON = 'SCC'
END_CODON   = 'SSS'
FUSION_MARKER = 'OCC'

POSITION_CODONS = [
    'COO', 'COS', 'CSC', 'CSN', 'CSO', 'CSS',
    'OOC', 'OON', 'OOO', 'OOS', 'SCO', 'SCS',
    'SNC', 'SNN', 'SNO', 'SNS',
]
RING_REF_CODONS = ['SOC', 'SON', 'SOO', 'SOS', 'SSC', 'SSN', 'SSO']

ATOM_PRIORITY = {
    'C': 0,
    'O': 1,
    'N': 2,
    'S': 3,
    'P': 4,
    'F': 5,
    'Cl': 6,
    'Br': 7,
    'I': 8,
    'B': 9,
}


def atom_priority(atom):
    return ATOM_PRIORITY.get(atom.GetSymbol(), 99)


def atom_order_key(mol, canon, atom_idx):
    atom = mol.GetAtomWithIdx(atom_idx)
    aromatic_rank = 0 if atom.GetIsAromatic() else 1
    return (atom_priority(atom), aromatic_rank, canon[atom_idx], atom_idx)


def neighbor_order_key(mol, canon, atom_idx, bond):
    atom = mol.GetAtomWithIdx(atom_idx)
    bond_order_rank = {
        Chem.BondType.AROMATIC: 0,
        Chem.BondType.DOUBLE: 1,
        Chem.BondType.TRIPLE: 2,
        Chem.BondType.SINGLE: 3,
    }.get(bond.GetBondType(), 9)
    return atom_order_key(mol, canon, atom_idx) + (bond_order_rank, bond.GetIdx())


def get_atom_annotations(atom, atom_h_counts=None):
    ann = []
    # Charge
    chg = atom.GetFormalCharge()
    if chg == 0: ann.append('CCX')
    elif chg == 1: ann.append('CXN')
    elif chg >= 2: ann.append('CXS')
    elif chg == -1: ann.append('CXO')
    elif chg <= -2: ann.append('CXX')

    # Stereo
    if atom.HasProp('_CIPCode'):
        cip = atom.GetProp('_CIPCode')
        if cip == 'R': ann.append('SXN')
        elif cip == 'S': ann.append('SXO')

    # Pharmacophore
    sym = atom.GetSymbol()
    if sym in ('N', 'O'):
        ann.append('OXN')  # HBA
        num_h = atom_h_counts[atom.GetIdx()] if atom_h_counts else atom.GetTotalNumHs()
        for _ in range(num_h):
            ann.append('OXO')  # Polar H (one per hydrogen)
    return ann

def _canonical_ez_label(bond, mol, canon):
    """
    Return 'E', 'Z', or None for a double bond, normalized so that the
    label is defined with respect to the highest-canonical-rank neighbor on
    each side. This matches the convention used by the decoder.
    """
    stereo = bond.GetStereo()
    if stereo not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
        return None
    sa = list(bond.GetStereoAtoms())
    if len(sa) != 2:
        return None
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    begin_nbrs = [n.GetIdx() for n in bond.GetBeginAtom().GetNeighbors() if n.GetIdx() != a2]
    end_nbrs = [n.GetIdx() for n in bond.GetEndAtom().GetNeighbors() if n.GetIdx() != a1]
    if not begin_nbrs or not end_nbrs:
        return None
    canonical_begin = max(begin_nbrs, key=lambda i: canon[i])
    canonical_end = max(end_nbrs, key=lambda i: canon[i])

    # Which side does each stereoAtom belong to?
    if sa[0] in begin_nbrs and sa[1] in end_nbrs:
        rdkit_begin, rdkit_end = sa[0], sa[1]
    elif sa[1] in begin_nbrs and sa[0] in end_nbrs:
        rdkit_begin, rdkit_end = sa[1], sa[0]
    else:
        return None  # Should not happen on a sane mol

    # Each side flip independently toggles E<->Z.
    flips = 0
    if rdkit_begin != canonical_begin:
        flips += 1
    if rdkit_end != canonical_end:
        flips += 1

    base = 'E' if stereo == Chem.BondStereo.STEREOE else 'Z'
    if flips % 2 == 1:
        base = 'Z' if base == 'E' else 'E'
    return base


def get_bond_annotations(bond, aromatic_bonds=None, mol=None, canon=None):
    ann = []
    is_aromatic = aromatic_bonds is not None and bond.GetIdx() in aromatic_bonds
    # Mobility
    if bond.IsInRing():
        if is_aromatic:
            ann.append('NXO')  # aromatic_locked
        else:
            ann.append('NXS')  # ring_constrained
    else:
        if bond.GetBondType() == Chem.BondType.SINGLE:
            b_a1 = bond.GetBeginAtom()
            b_a2 = bond.GetEndAtom()
            if b_a1.GetDegree() == 1 or b_a2.GetDegree() == 1:
                ann.append('NCX')  # non_rotatable if terminal
            else:
                ann.append('NXC')  # rotatable
        else:
            ann.append('NCX')  # double/triple non rotatable

    # Stereo (normalized to canonical-rank-max neighbor convention)
    if mol is not None and canon is not None:
        ez = _canonical_ez_label(bond, mol, canon)
        if ez == 'E':
            ann.append('SOX')
        elif ez == 'Z':
            ann.append('SNX')
    else:
        # Fallback (raw RDKit stereo) — kept for backward compatibility,
        # but encode() now always passes mol+canon so we go through the
        # canonical branch above.
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOE:
            ann.append('SOX')
        elif stereo == Chem.BondStereo.STEREOZ:
            ann.append('SNX')

    return ann


def encode(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise EncodeError(f"Invalid SMILES: {smiles}")
    if mol.GetNumAtoms() == 0:
        raise EncodeError("Empty molecule is not encodable")
    atom_h_counts = {a.GetIdx(): a.GetTotalNumHs() for a in mol.GetAtoms()}
    aromatic_bonds = {b.GetIdx() for b in mol.GetBonds() if b.GetIsAromatic()}
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        raise EncodeError(f"Kekulization failed: {e}") from e
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    canon = list(Chem.CanonicalRankAtoms(mol))
    all_rings = mol.GetRingInfo().AtomRings()
    a2r = defaultdict(list)
    for ri, ring in enumerate(all_rings):
        for ai in ring:
            a2r[ai].append(ri)

    adj = defaultdict(list)
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj[a1].append((a2, b)); adj[a2].append((a1, b))

    codons = [START_CODON]
    vis_a, vis_r = set(), set()
    _slot_next = [0]
    _slot_free = []
    def alloc():
        if _slot_free:
            return _slot_free.pop(0)
        s = _slot_next[0]
        if s >= len(RING_OPEN):
            raise EncodeError(
                f"Ring label capacity exceeded: need slot {s}, max is {len(RING_OPEN) - 1}"
            )
        _slot_next[0] += 1
        return s
    def free(s):
        if s not in _slot_free:
            _slot_free.append(s)

    completed_rings = []
    def uring(ai):
        for ri in a2r.get(ai, []):
            if ri not in vis_r: return ri
        return None

    def neighbor_key(cur, nb, bond):
        return neighbor_order_key(mol, canon, nb, bond)

    terminals = [a.GetIdx() for a in mol.GetAtoms() if a.GetDegree() == 1]
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
        if not bond:
            return
        bond_type = bond.GetBondType()
        if bond_type not in BOND_TYPE_TO_CODON:
            raise EncodeError(f"Unsupported bond type: {bond_type}")
        codons.append(BOND_TYPE_TO_CODON[bond_type])
        codons.extend(get_bond_annotations(bond, aromatic_bonds, mol=mol, canon=canon))

    def emit_atom(ai):
        atom = mol.GetAtomWithIdx(ai)
        symbol = atom.GetSymbol()
        if symbol not in ATOM_TO_CODON:
            raise EncodeError(f"Unsupported atom type: {symbol}")
        codons.append(ATOM_TO_CODON[symbol])
        codons.extend(get_atom_annotations(atom, atom_h_counts))

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
                    codons.append(BRANCH_OPEN[d])
                    traverse(n, b)
                    codons.append(BRANCH_CLOSE[d])
                else:
                    traverse(n, b)

    def enc_ring(ri, entry, from_bond):
        if ri in vis_r: return
        vis_r.add(ri)
        slot = alloc()
        if from_bond is not None: emit_bond(from_bond)
        codons.append(RING_OPEN[slot])

        ring_atoms = rorder(ri, entry)
        ring_set = set(ring_atoms)
        ring_atom_list = []
        deferred = []
        branch_label_counter = 0
        seen_fused_targets = set()
        
        prev = None
        for i, a in enumerate(ring_atoms):
            already_visited = a in vis_a
            vis_a.add(a)
            if i > 0 and prev is not None:
                b = mol.GetBondBetweenAtoms(prev, a)
                if b: emit_bond(b)
            
            if already_visited:
                codons.append(FUSION_MARKER)
                fusion_written = False
                for cri, cring in enumerate(completed_rings):
                    if a in cring:
                        pos = cring.index(a)
                        if cri >= len(RING_REF_CODONS):
                            raise EncodeError(
                                f"Completed ring reference capacity exceeded for fused atom {a}: ring index {cri}"
                            )
                        if pos >= len(POSITION_CODONS):
                            raise EncodeError(
                                f"Position codon capacity exceeded for fused atom {a}: position {pos}"
                            )
                        codons.append(RING_REF_CODONS[cri])
                        codons.append(POSITION_CODONS[pos])
                        fusion_written = True
                        break
                if not fusion_written:
                    raise EncodeError(f"Failed to resolve fused atom reference for atom {a}")
                # Fused atom: only symbol, no annotations (already written on first visit)
                atom = mol.GetAtomWithIdx(a)
                symbol = atom.GetSymbol()
                if symbol not in ATOM_TO_CODON:
                    raise EncodeError(f"Unsupported fused atom type: {symbol}")
                codons.append(ATOM_TO_CODON[symbol])
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
                            codons.append(RING_OPEN[fused_label_slot])
                            deferred.append(('fused', fused_label_slot, a, target_ring, nb, nb_bond))
                    else:
                        br_idx = branch_label_counter % len(BRANCH_OPEN)
                        codons.append(BRANCH_OPEN[br_idx])
                        deferred.append(('branch', br_idx, a, nb, nb_bond))
                        branch_label_counter += 1
            prev = a

        closing_bond = mol.GetBondBetweenAtoms(ring_atoms[-1], ring_atoms[0])
        if closing_bond: emit_bond(closing_bond)
        codons.append(RING_CLOSE[slot])
        free(slot)
        completed_rings.append(ring_atom_list)
        
        for item in deferred:
            if item[0] == 'branch':
                _, br_idx, from_atom, nb, nb_bond = item
                if nb not in vis_a:
                    codons.append(BRANCH_OPEN[br_idx])
                    traverse(nb, nb_bond)
                    codons.append(BRANCH_CLOSE[br_idx])
                else:
                    codons.append(BRANCH_OPEN[br_idx])
                    codons.append(BRANCH_CLOSE[br_idx])
            elif item[0] == 'fused':
                _, fused_label_slot, from_atom, target_ring, nb, nb_bond = item
                if target_ring not in vis_r:
                    parent_ring_idx = None
                    pos = None
                    for cri, cring in enumerate(completed_rings):
                        if from_atom in cring:
                            parent_ring_idx = cri
                            pos = cring.index(from_atom)
                            break
                    if parent_ring_idx is None or pos is None:
                        raise EncodeError(
                            f"Failed to resolve deferred fused parent reference for atom {from_atom}"
                        )
                    if parent_ring_idx >= len(RING_REF_CODONS):
                        raise EncodeError(
                            f"Deferred fused parent ring reference capacity exceeded: ring index {parent_ring_idx}"
                        )
                    if pos >= len(POSITION_CODONS):
                        raise EncodeError(
                            f"Deferred fused parent position capacity exceeded: position {pos}"
                        )
                    codons.append(RING_OPEN[fused_label_slot])
                    codons.append(RING_REF_CODONS[parent_ring_idx])
                    codons.append(POSITION_CODONS[pos])
                    _slot_free.insert(0, fused_label_slot)
                    actual_bond = mol.GetBondBetweenAtoms(from_atom, nb)
                    if actual_bond is None:
                        raise EncodeError(
                            f"Missing bond for deferred fused transition {from_atom}->{nb}"
                        )
                    enc_ring(target_ring, nb, actual_bond)
                else:
                    free(fused_label_slot)

    traverse(start)
    codons.append(END_CODON)
    return codons

def main():
    parser = argparse.ArgumentParser(description="Encode a SMILES string into a MolCodon sequence.")
    parser.add_argument("smiles", help="Input SMILES string")
    parser.add_argument(
        "--one-per-line",
        action="store_true",
        help="Print one codon per line instead of a whitespace-separated sequence",
    )
    args = parser.parse_args()

    codons = encode(args.smiles)
    if args.one_per_line:
        print("\n".join(codons))
    else:
        print(" ".join(codons))

if __name__ == "__main__":
    main()
