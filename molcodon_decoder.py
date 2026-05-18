#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem

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
CODON_TO_BOND = {v: k for k, v in BOND_TYPE_TO_CODON.items()}

BRANCH_OPEN = ['NNC', 'NOC', 'NOS', 'NSC']
BRANCH_CLOSE = ['NNN', 'NON', 'NOO', 'NSN']
BRANCH_CLOSE_BY_OPEN = dict(zip(BRANCH_OPEN, BRANCH_CLOSE))

RING_OPEN = ['NNO', 'NSO', 'OSO', 'OCN', 'OCO', 'OCS', 'ONC', 'ONN']
RING_CLOSE = ['NNS', 'NSS', 'OSS', 'ONO', 'ONS', 'OSC', 'OSN', 'SCN']
RING_CLOSE_BY_OPEN = dict(zip(RING_OPEN, RING_CLOSE))
RING_OPEN_BY_CLOSE = dict(zip(RING_CLOSE, RING_OPEN))

START_CODON = 'SCC'
END_CODON = 'SSS'
FUSION_MARKER = 'OCC'

POSITION_CODONS = [
    'COO', 'COS', 'CSC', 'CSN', 'CSO', 'CSS',
    'OOC', 'OON', 'OOO', 'OOS', 'SCO', 'SCS',
    'SNC', 'SNN', 'SNO', 'SNS',
]
RING_REF_CODONS = ['SOC', 'SON', 'SOO', 'SOS', 'SSC', 'SSN', 'SSO']
POSITION_BY_CODON = {c: i for i, c in enumerate(POSITION_CODONS)}
RING_REF_BY_CODON = {c: i for i, c in enumerate(RING_REF_CODONS)}

ATOM_ANNOTATIONS = {
    'CCX': ('charge', 0), 'CXN': ('charge', 1), 'CXS': ('charge', 2), 'CXO': ('charge', -1), 'CXX': ('charge', -2),
    'SXN': ('atom_stereo', 'R'), 'SXO': ('atom_stereo', 'S'),
    'OXN': ('pharmacophore', 'hbond_acceptor'), 'OXO': ('pharmacophore', 'polar_h'),
}
BOND_ANNOTATIONS = {
    'NXC': ('bond_mobility', 'rotatable'), 'NCX': ('bond_mobility', 'non_rotatable'), 'NXS': ('bond_mobility', 'ring_constrained'),
    'NXO': ('bond_mobility', 'aromatic_locked'), 'SOX': ('bond_stereo', 'E'), 'SNX': ('bond_stereo', 'Z'),
}


class DecodeError(Exception):
    pass


@dataclass
class PendingBond:
    bond_type: Chem.BondType
    annotations: Dict[str, object] = field(default_factory=dict)


class MolcodonV2Decoder:
    def __init__(self, tokens: List[str]):
        self.tokens = tokens
        self.pos = 0
        self.mol = Chem.RWMol()
        self.completed_rings: List[List[int]] = []
        self.atom_props: Dict[int, Dict[str, object]] = {}
        self.bond_props: Dict[Tuple[int, int], Dict[str, object]] = {}

    def peek(self) -> Optional[str]:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def pop(self) -> str:
        tok = self.peek()
        if tok is None:
            raise DecodeError('Unexpected end of token stream')
        self.pos += 1
        return tok

    def expect(self, token: str) -> None:
        got = self.pop()
        if got != token:
            raise DecodeError(f'Expected {token}, got {got} at token {self.pos - 1}')

    def parse_atom_annotations(self) -> Dict[str, object]:
        anns: Dict[str, object] = {'pharmacophore': []}
        while self.peek() in ATOM_ANNOTATIONS:
            kind, value = ATOM_ANNOTATIONS[self.pop()]
            if kind == 'pharmacophore':
                anns['pharmacophore'].append(value)
            else:
                anns[kind] = value
        return anns

    def parse_bond(self) -> PendingBond:
        tok = self.pop()
        if tok not in CODON_TO_BOND:
            raise DecodeError(f'Expected bond codon, got {tok} at token {self.pos - 1}')
        anns: Dict[str, object] = {}
        while self.peek() in BOND_ANNOTATIONS:
            kind, value = BOND_ANNOTATIONS[self.pop()]
            anns[kind] = value
        return PendingBond(CODON_TO_BOND[tok], anns)

    def add_atom(self, symbol: str, annotations: Dict[str, object]) -> int:
        atom = Chem.Atom(symbol)
        charge = annotations.get('charge')
        if isinstance(charge, int):
            atom.SetFormalCharge(charge)
        idx = self.mol.AddAtom(atom)
        self.atom_props[idx] = annotations
        return idx

    def set_bond_props(self, a1: int, a2: int, annotations: Dict[str, object]) -> None:
        key = tuple(sorted((a1, a2)))
        self.bond_props[key] = annotations

    def add_bond(self, a1: int, a2: int, pending: PendingBond) -> None:
        if a1 == a2:
            return
        existing = self.mol.GetBondBetweenAtoms(a1, a2)
        if existing is None:
            self.mol.AddBond(a1, a2, pending.bond_type)
        self.set_bond_props(a1, a2, pending.annotations)

    def parse_atom_or_fusion(self) -> Tuple[int, bool]:
        if self.peek() == FUSION_MARKER:
            self.pop()
            ring_ref_tok = self.pop()
            pos_ref_tok = self.pop()
            ring_idx = RING_REF_BY_CODON.get(ring_ref_tok)
            pos_idx = POSITION_BY_CODON.get(pos_ref_tok)
            if ring_idx is None or pos_idx is None:
                raise DecodeError('Invalid fusion reference')
            try:
                atom_idx = self.completed_rings[ring_idx][pos_idx]
            except IndexError as e:
                raise DecodeError('Fusion reference points outside completed ring table') from e
            atom_codon = self.pop()
            if atom_codon not in CODON_TO_ATOM:
                raise DecodeError(f'Expected atom codon after fusion marker, got {atom_codon}')
            if self.peek() in ATOM_ANNOTATIONS:
                raise DecodeError('Fusion atom must not carry atom annotations')
            return atom_idx, True
        atom_codon = self.pop()
        if atom_codon not in CODON_TO_ATOM:
            raise DecodeError(f'Expected atom codon, got {atom_codon} at token {self.pos - 1}')
        atom_idx = self.add_atom(CODON_TO_ATOM[atom_codon], self.parse_atom_annotations())
        return atom_idx, False

    def parse_branch_expansion(self, root_atom: int, open_tok: str) -> None:
        self.expect(open_tok)
        close_tok = BRANCH_CLOSE_BY_OPEN[open_tok]
        if self.peek() == close_tok:
            self.pop()
            return
        pending = self.parse_bond()
        child = self.parse_node(root_atom, pending)
        self.parse_followups(child)
        self.expect(close_tok)

    def parse_ring(self, parent_atom, incoming, ring_open_tok: str) -> int:
        self.expect(ring_open_tok)
        first_atom, fused = self.parse_atom_or_fusion()
        if incoming is not None and parent_atom is not None:
            self.add_bond(parent_atom, first_atom, incoming)
        ring_atoms = [first_atom]
        current = first_atom
        deferred: List[Tuple[str, str, int]] = []

        close_tok = RING_CLOSE_BY_OPEN[ring_open_tok]
        while True:
            tok = self.peek()
            if tok is None:
                raise DecodeError('Unterminated ring')
            if tok in BRANCH_OPEN:
                deferred.append(('branch', self.pop(), current))
                continue
            if tok in RING_OPEN:
                deferred.append(('fused', self.pop(), current))
                continue
            if tok == close_tok:
                self.pop()
                break
            pending = self.parse_bond()
            if self.peek() == close_tok:
                self.add_bond(current, first_atom, pending)
                self.pop()
                break
            nxt, fused = self.parse_atom_or_fusion()
            self.add_bond(current, nxt, pending)
            current = nxt
            ring_atoms.append(nxt)

        self.completed_rings.append(ring_atoms)

        for kind, marker, root in deferred:
            if kind == 'branch':
                self.parse_branch_expansion(root, marker)
            else:
                open_tok = self.peek()
                # A deferred fused placeholder may become redundant if that target
                # ring was already visited indirectly while replaying an earlier
                # deferred fused expansion. In that case the encoder emits no replay
                # block for this marker, so the decoder must skip it as well.
                if open_tok != marker:
                    continue
                self.pop()
                ring_ref_tok = self.pop()
                pos_ref_tok = self.pop()
                if ring_ref_tok not in RING_REF_BY_CODON:
                    raise DecodeError('Deferred fused ring is missing ring reference token')
                if pos_ref_tok not in POSITION_BY_CODON:
                    raise DecodeError('Deferred fused ring is missing position token')
                bond = self.parse_bond()
                self.parse_ring(root, bond, marker)
        return current

    def parse_node(self, parent_atom: int, incoming: PendingBond) -> int:
        tok = self.peek()
        if tok in RING_OPEN:
            return self.parse_ring(parent_atom, incoming, tok)
        atom_idx, _ = self.parse_atom_or_fusion()
        self.add_bond(parent_atom, atom_idx, incoming)
        return atom_idx

    def parse_followups(self, current_atom: int) -> int:
        while True:
            tok = self.peek()
            if tok in BRANCH_OPEN:
                self.parse_branch_expansion(current_atom, tok)
                continue
            if tok in CODON_TO_BOND:
                pending = self.parse_bond()
                current_atom = self.parse_node(current_atom, pending)
                continue
            break
        return current_atom

    def finalize(self) -> Chem.Mol:
        mol = self.mol.GetMol()
        # Store annotations as custom properties for downstream consumers AND
        # collect the targets we will use to actually set RDKit-level stereo.
        target_atom_stereo: Dict[int, str] = {}   # atom_idx -> 'R'/'S'
        target_bond_stereo: Dict[Tuple[int, int], str] = {}  # (a1,a2) -> 'E'/'Z'

        for idx, anns in self.atom_props.items():
            atom = mol.GetAtomWithIdx(idx)
            if anns.get('atom_stereo'):
                atom.SetProp('molcodon_stereo', anns['atom_stereo'])
                target_atom_stereo[idx] = anns['atom_stereo']
            if anns.get('pharmacophore'):
                atom.SetProp('molcodon_pharmacophore', ','.join(anns['pharmacophore']))

        for (a1, a2), anns in self.bond_props.items():
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None:
                continue
            if anns.get('bond_mobility'):
                bond.SetProp('molcodon_bond_mobility', str(anns['bond_mobility']))
            if anns.get('bond_stereo') in ('E', 'Z'):
                target_bond_stereo[(a1, a2)] = anns['bond_stereo']

        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            raise DecodeError(f'Sanitize failed: {e}') from e

        # --- Restore atom (tetrahedral) stereo via R/S target matching ---
        # ChiralTag (CW/CCW) is neighbor-order-dependent, so we cannot map
        # R/S -> ChiralTag directly. Strategy: try CW, ask RDKit what CIP
        # it produces, flip to CCW if it doesn't match the target.
        if target_atom_stereo:
            for atom_idx, target_rs in target_atom_stereo.items():
                atom = mol.GetAtomWithIdx(atom_idx)
                # Try CW first
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                produced = atom.GetPropsAsDict().get('_CIPCode', None)
                if produced != target_rs:
                    atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                    produced = atom.GetPropsAsDict().get('_CIPCode', None)
                    if produced != target_rs:
                        # Either the atom has no real stereo (e.g. degenerate
                        # neighbors after sanitize), or RDKit refuses to assign.
                        # Clear the tag rather than leave bogus chirality.
                        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        # --- Restore bond (double-bond) stereo via E/Z target matching ---
        # E/Z depends on which two neighbors are chosen as stereo atoms.
        # We deterministically pick the highest-canonical-rank neighbor on each
        # side. The encoder uses the SAME rule when emitting SOX/SNX, so the
        # E/Z label is invariant under this convention.
        #
        # Implementation note: we set StereoAtoms + STEREOE/Z and then call
        # AssignStereochemistry with cleanIt=False. cleanIt=True would re-derive
        # stereo from BondDir flags (which we don't set) and wipe it. Downstream
        # consumers that call AssignStereochemistry(cleanIt=True) on this mol
        # will lose the stereo unless they preserve it themselves; the SMILES
        # writer reads stereo directly so MolToSmiles(isomericSmiles=True)
        # still produces the correct directional SMILES.
        if target_bond_stereo:
            canon = list(Chem.CanonicalRankAtoms(mol, breakTies=True))
            for (a1, a2), target_ez in target_bond_stereo.items():
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                begin_atom = mol.GetAtomWithIdx(a1)
                end_atom = mol.GetAtomWithIdx(a2)
                begin_nbrs = [n.GetIdx() for n in begin_atom.GetNeighbors() if n.GetIdx() != a2]
                end_nbrs = [n.GetIdx() for n in end_atom.GetNeighbors() if n.GetIdx() != a1]
                if not begin_nbrs or not end_nbrs:
                    continue
                begin_ref = max(begin_nbrs, key=lambda i: canon[i])
                end_ref = max(end_nbrs, key=lambda i: canon[i])
                bond.SetStereoAtoms(begin_ref, end_ref)
                if target_ez == 'E':
                    bond.SetStereo(Chem.BondStereo.STEREOE)
                else:
                    bond.SetStereo(Chem.BondStereo.STEREOZ)
            # cleanIt=False so we don't lose the stereo we just set.
            Chem.AssignStereochemistry(mol, cleanIt=False, force=True)

        return mol

    def decode(self) -> Chem.Mol:
        self.expect(START_CODON)
        tok = self.peek()
        if tok in RING_OPEN:
            # parse_ring calls self.expect(ring_open_tok) internally — don't pop here
            start_idx = self.parse_ring(None, None, tok)
        elif tok in CODON_TO_ATOM:
            self.pop()
            start_idx = self.add_atom(CODON_TO_ATOM[tok], self.parse_atom_annotations())
        else:
            raise DecodeError(f'Expected atom or ring open after SCC, got {tok}')
        self.parse_followups(start_idx)
        self.expect(END_CODON)
        if self.pos != len(self.tokens):
            raise DecodeError('Trailing tokens remain after END codon')
        return self.finalize()


def decode_sequence(sequence: str) -> Chem.Mol:
    tokens = [tok for tok in sequence.split() if tok]
    return MolcodonV2Decoder(tokens).decode()


def canonicalize_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smiles}')
    return Chem.MolToSmiles(mol, canonical=True)


def main() -> None:
    parser = argparse.ArgumentParser(description='Decode MolCodon v2 sequence back to an RDKit molecule.')
    parser.add_argument('input', help='Whitespace-separated MolCodon sequence or a file containing one')
    parser.add_argument('--from-file', action='store_true', help='Treat input as path to file containing the sequence')
    args = parser.parse_args()

    sequence = open(args.input).read().strip() if args.from_file else args.input
    mol = decode_sequence(sequence)
    print(Chem.MolToSmiles(mol, canonical=True))


if __name__ == '__main__':
    main()
