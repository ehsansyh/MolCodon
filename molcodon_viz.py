#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from molcodon_match import MatchResult, match_smiles

CODON_CLASS = {
    'SCC': 'start', 'SSS': 'end',
    'NCC': 'bond', 'NCN': 'bond', 'NCO': 'bond', 'NCS': 'bond',
    'CCC': 'atom', 'CCN': 'atom', 'CCO': 'atom', 'CCS': 'atom', 'CNC': 'atom', 'CNN': 'atom', 'CNO': 'atom', 'CNS': 'atom', 'COC': 'atom', 'CON': 'atom',
    'NNC': 'branch-open', 'NOC': 'branch-open', 'NOS': 'branch-open', 'NSC': 'branch-open',
    'NNN': 'branch-close', 'NON': 'branch-close', 'NOO': 'branch-close', 'NSN': 'branch-close',
    'NNO': 'ring-open', 'NSO': 'ring-open', 'OSO': 'ring-open', 'OCN': 'ring-open', 'OCO': 'ring-open', 'OCS': 'ring-open', 'ONC': 'ring-open', 'ONN': 'ring-open',
    'NNS': 'ring-close', 'NSS': 'ring-close', 'OSS': 'ring-close', 'ONO': 'ring-close', 'ONS': 'ring-close', 'OSC': 'ring-close', 'OSN': 'ring-close', 'SCN': 'ring-close',
    'OXN': 'pharm', 'OXO': 'pharm', 'CXN': 'pharm', 'CXS': 'pharm', 'CXO': 'pharm', 'CXX': 'pharm',
    'NXO': 'bond-ann', 'NXS': 'bond-ann', 'NXC': 'bond-ann', 'NCX': 'bond-ann', 'SOX': 'stereo', 'SNX': 'stereo', 'SXN': 'stereo', 'SXO': 'stereo',
    'OCC': 'fusion',
    'SOC': 'ring-ref', 'SON': 'ring-ref', 'SOO': 'ring-ref', 'SOS': 'ring-ref', 'SSC': 'ring-ref', 'SSN': 'ring-ref', 'SSO': 'ring-ref',
    'COO': 'position', 'COS': 'position', 'CSC': 'position', 'CSN': 'position', 'CSO': 'position', 'CSS': 'position', 'OOC': 'position', 'OON': 'position', 'OOO': 'position', 'OOS': 'position', 'SCO': 'position', 'SCS': 'position', 'SNC': 'position', 'SNN': 'position', 'SNO': 'position', 'SNS': 'position',
}

# ── Colour palette ──────────────────────────────────────────────────
# Each component type gets a base colour.
# Pairs within the same type get the SAME colour so identical
# structural elements are visually linked across ref and hit.
# Colours are deliberately light/pastel so they layer readably.

MATCH_PALETTE = {
    'ring':       (0.40, 0.78, 0.94),   # soft cyan
    'branch':     (0.96, 0.68, 0.38),   # soft orange
    'attachment': (0.74, 0.56, 0.92),   # soft violet
    'pharm':      (0.48, 0.86, 0.52),   # soft green
}

# CSS hex equivalents (for the HTML tables / swatches)
MATCH_PALETTE_HEX = {
    'ring':       '#66c7f0',
    'branch':     '#f5ad61',
    'attachment': '#bd8feb',
    'pharm':      '#7adb85',
}


def _add_opacity_to_svg(svg: str, opacity: float = 0.45) -> str:
    """Post-process RDKit SVG to make highlight fills semi-transparent."""
    svg = re.sub(
        r'fill-opacity:1(?!\.)',
        f'fill-opacity:{opacity}',
        svg,
    )
    return svg


def _build_highlight_maps(matched_components: dict, side: str) -> Tuple[Dict[int, Tuple], Dict[int, Tuple]]:
    """Build atom→colour and bond→colour dicts for one side (ref / hit).

    Priority (last wins): ring < branch < pharm
    so pharmacophore annotations paint over structural highlights.
    """
    atoms: Dict[int, Tuple] = {}
    bonds: Dict[int, Tuple] = {}

    def _add(alist, blist, colour):
        for a in alist:
            atoms[a] = colour
        for b in blist:
            bonds[b] = colour

    key_a = 'ref_atoms' if side == 'ref' else 'hit_atoms'
    key_b = 'ref_bonds' if side == 'ref' else 'hit_bonds'

    c = MATCH_PALETTE['ring']
    for r in matched_components.get('rings', []):
        _add(r.get(key_a, []), r.get(key_b, []), c)

    c = MATCH_PALETTE['branch']
    for b in matched_components.get('branches', []):
        _add(b.get(key_a, []), b.get(key_b, []), c)

    c = MATCH_PALETTE['pharm']
    key_atom = 'ref_atom' if side == 'ref' else 'hit_atom'
    for p in matched_components.get('pharmacophores', []):
        a = p.get(key_atom)
        if a is not None:
            atoms[a] = c

    c = MATCH_PALETTE['attachment']
    key_att = 'ref_attach_atom' if side == 'ref' else 'hit_attach_atom'
    for att in matched_components.get('attachments', []):
        a = att.get(key_att)
        if a is not None:
            atoms[a] = c

    return atoms, bonds


def mol_to_svg(mol: Chem.Mol, highlight_atoms: Dict[int, Tuple],
               highlight_bonds: Dict[int, Tuple],
               width: int = 420, height: int = 320) -> str:
    mol = Chem.Mol(mol)
    rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addAtomIndices = True
    opts.highlightBondWidthMultiplier = 16
    drawer.DrawMolecule(
        mol,
        highlightAtoms=list(highlight_atoms.keys()),
        highlightAtomColors=highlight_atoms,
        highlightBonds=list(highlight_bonds.keys()),
        highlightBondColors=highlight_bonds,
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return _add_opacity_to_svg(svg)


def render_token_sequence(tokens: List[str], spans: Iterable[Tuple[int, int]]) -> str:
    matched = set()
    for start, end in spans:
        matched.update(range(start, end + 1))
    pieces = []
    for idx, tok in enumerate(tokens):
        css = CODON_CLASS.get(tok, 'other')
        classes = ['codon', css]
        if idx in matched:
            classes.append('matched')
        pieces.append(f'<span class="{" ".join(classes)}" title="token {idx}">{html.escape(tok)}</span>')
    return ' '.join(pieces)


def component_table(scores: dict, tanimoto: float) -> str:
    rows = [
        ('Overall', f"{scores['overall']:.2f}", ''),
        ('Ring F1', f"{scores['ring_f1']:.2f} ({scores['ring_match']})", MATCH_PALETTE_HEX['ring']),
        ('Branch F1', f"{scores['branch_f1']:.2f} ({scores['branch_match']})", MATCH_PALETTE_HEX['branch']),
        ('Attachment F1', f"{scores['attachment_f1']:.2f} ({scores['attachment_match']})", MATCH_PALETTE_HEX['attachment']),
        ('Pharmacophore F1', f"{scores['pharmacophore_f1']:.2f} ({scores['pharmacophore_match']})", MATCH_PALETTE_HEX['pharm']),
        ('Tanimoto (ref)', f"{tanimoto:.2f}", ''),
    ]
    body = ''
    for label, value, colour in rows:
        swatch = f'<span style="display:inline-block;width:14px;height:14px;border-radius:4px;background:{colour};margin-right:8px;vertical-align:middle;"></span>' if colour else ''
        body += f'<tr><th>{swatch}{html.escape(label)}</th><td>{html.escape(value)}</td></tr>'
    return f'<table class="scorecard">{body}</table>'


def legend_html() -> str:
    items = [
        (MATCH_PALETTE_HEX['ring'], 'Ring match'),
        (MATCH_PALETTE_HEX['branch'], 'Branch match'),
        (MATCH_PALETTE_HEX['attachment'], 'Attachment match'),
        (MATCH_PALETTE_HEX['pharm'], 'Pharmacophore match'),
    ]
    out = []
    for colour, label in items:
        out.append(f'<div class="legend-item"><span class="legend-swatch" style="background:{colour};"></span>{html.escape(label)}</div>')
    return '<div class="legend">' + ''.join(out) + '</div>'


def _atoms_str(atoms) -> str:
    if not atoms:
        return '-'
    return ', '.join(str(a) for a in atoms)


def build_components_tables(matched_components: dict) -> str:
    html_out = []

    # ── Rings ────────────────────────────────────────────
    rings = matched_components.get('rings', [])
    if rings:
        c = MATCH_PALETTE_HEX['ring']
        html_out.append(f'<h3><span class="swatch" style="background:{c}"></span> Matched Rings</h3>')
        html_out.append('<table class="comp-table"><tr><th>Ref Ring</th><th>Hit Ring</th><th>Ref Atoms</th><th>Hit Atoms</th></tr>')
        for r in rings:
            html_out.append(
                f'<tr><td>Ring {r["ref_id"]}</td><td>Ring {r["hit_id"]}</td>'
                f'<td>{_atoms_str(r.get("ref_atoms"))}</td>'
                f'<td>{_atoms_str(r.get("hit_atoms"))}</td></tr>'
            )
        html_out.append('</table>')

    # ── Branches ─────────────────────────────────────────
    branches = matched_components.get('branches', [])
    if branches:
        c = MATCH_PALETTE_HEX['branch']
        html_out.append(f'<h3><span class="swatch" style="background:{c}"></span> Matched Branches</h3>')
        html_out.append('<table class="comp-table"><tr><th>Ref Branch</th><th>Hit Branch</th><th>Parent</th><th>Ref Atoms</th><th>Hit Atoms</th></tr>')
        for b in branches:
            html_out.append(
                f'<tr><td>Branch {b["ref_id"]}</td><td>Branch {b["hit_id"]}</td>'
                f'<td>{html.escape(str(b.get("parent_type", "")))}</td>'
                f'<td>{_atoms_str(b.get("ref_atoms"))}</td>'
                f'<td>{_atoms_str(b.get("hit_atoms"))}</td></tr>'
            )
        html_out.append('</table>')

    # ── Attachments ──────────────────────────────────────
    attachments = matched_components.get('attachments', [])
    if attachments:
        c = MATCH_PALETTE_HEX['attachment']
        html_out.append(f'<h3><span class="swatch" style="background:{c}"></span> Matched Attachments</h3>')
        html_out.append('<table class="comp-table"><tr><th>Ref Branch</th><th>Hit Branch</th><th>Parent Type</th><th>Ref Atom</th><th>Hit Atom</th><th>Similarity</th></tr>')
        for a in attachments:
            sim = a.get("similarity", 0)
            sim_pct = f'{sim * 100:.0f}%'
            html_out.append(
                f'<tr><td>Branch {a["ref_branch_id"]}</td><td>Branch {a["hit_branch_id"]}</td>'
                f'<td>{html.escape(a["parent_type"])}</td>'
                f'<td>{a.get("ref_attach_atom", "-")}</td>'
                f'<td>{a.get("hit_attach_atom", "-")}</td>'
                f'<td>{html.escape(sim_pct)}</td></tr>'
            )
        html_out.append('</table>')

    # ── Pharmacophores ───────────────────────────────────
    pharms = matched_components.get('pharmacophores', [])
    if pharms:
        c = MATCH_PALETTE_HEX['pharm']
        html_out.append(f'<h3><span class="swatch" style="background:{c}"></span> Matched Pharmacophores</h3>')
        html_out.append('<table class="comp-table"><tr><th>Feature</th><th>Context</th><th>Ref Atom</th><th>Hit Atom</th></tr>')
        for p in pharms:
            html_out.append(
                f'<tr><td>{html.escape(p["label"])}</td><td>{html.escape(p["context"])}</td>'
                f'<td>{p["ref_atom"]}</td><td>{p["hit_atom"]}</td></tr>'
            )
        html_out.append('</table>')

    if not html_out:
        return '<p>No mapped components found.</p>'
    return ''.join(html_out)


def build_report(match: MatchResult, output_path: Path,
                 title: str = 'MolCodon match report') -> None:

    ref_atoms, ref_bonds = _build_highlight_maps(match.matched_components, 'ref')
    hit_atoms, hit_bonds = _build_highlight_maps(match.matched_components, 'hit')

    ref_svg = mol_to_svg(match.reference.mol, ref_atoms, ref_bonds)
    hit_svg = mol_to_svg(match.hit.mol, hit_atoms, hit_bonds)
    ref_tokens = render_token_sequence(match.reference.tokens, match.matched_ref_token_spans)
    hit_tokens = render_token_sequence(match.hit.tokens, match.matched_hit_token_spans)
    scorecard = component_table(match.scores, match.tanimoto)
    components_html = build_components_tables(match.matched_components)

    html_doc = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>{html.escape(title)}</title>
<style>
body {{ font-family: Inter, Arial, sans-serif; margin: 24px; color: #0f172a; background: #f8fafc; }}
.grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; align-items: start; }}
.panel {{ background: white; border-radius: 16px; padding: 16px; box-shadow: 0 8px 24px rgba(15,23,42,0.08); }}
h1, h2, h3 {{ margin: 0 0 12px 0; }}
h3 {{ display: flex; align-items: center; gap: 8px; }}
.swatch {{ display: inline-block; width: 16px; height: 16px; border-radius: 4px; flex-shrink: 0; }}
.structure svg {{ width: 100%; height: auto; border: 1px solid #e2e8f0; border-radius: 12px; background: #fff; }}
.tokens {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; line-height: 2.0; word-break: break-word; }}
.codon {{ display: inline-block; padding: 2px 6px; margin: 2px; border-radius: 999px; color: white; font-size: 13px; }}
.codon.matched {{ box-shadow: 0 0 0 3px rgba(34,197,94,0.25); transform: translateY(-1px); }}
.scorecard {{ border-collapse: collapse; width: 100%; }}
.scorecard th, .scorecard td {{ border-bottom: 1px solid #e2e8f0; padding: 10px 12px; text-align: left; }}
.scorecard th {{ width: 240px; color: #334155; }}
.comp-table {{ border-collapse: collapse; width: 100%; margin-bottom: 16px; }}
.comp-table th, .comp-table td {{ border-bottom: 1px solid #e2e8f0; padding: 8px 12px; text-align: left; font-size: 14px; }}
.comp-table th {{ color: #334155; font-weight: 600; background: #f8fafc; }}
.legend {{ display: flex; gap: 14px; flex-wrap: wrap; margin-top: 12px; }}
.legend-item {{ display: flex; align-items: center; gap: 8px; color: #334155; font-size: 14px; }}
.legend-swatch {{ width: 14px; height: 14px; border-radius: 4px; display: inline-block; opacity: 0.7; }}
.start, .end {{ background: #111827; }}
.atom {{ background: #1d4ed8; }}
.bond {{ background: #7c3aed; }}
.branch-open {{ background: #ea580c; }}
.branch-close {{ background: #c2410c; }}
.ring-open {{ background: #0891b2; }}
.ring-close {{ background: #0e7490; }}
.pharm {{ background: #16a34a; }}
.bond-ann {{ background: #9333ea; }}
.stereo {{ background: #be123c; }}
.fusion {{ background: #475569; }}
.ring-ref {{ background: #b45309; }}
.position {{ background: #64748b; }}
.other {{ background: #334155; }}
</style>
</head>
<body>
  <div class="panel">
    <h1>{html.escape(title)}</h1>
    <p>Matched structural components are highlighted with transparent colours.
       Atom indices on the structures correspond to the atom IDs in the tables below.</p>
    {scorecard}
    {legend_html()}
  </div>
  <div class="grid" style="margin-top:20px;">
    <div class="panel">
      <h2>Reference</h2>
      <div><strong>SMILES:</strong> {html.escape(match.reference.smiles)}</div>
      <div class="structure">{ref_svg}</div>
      <h3>MolCodon sequence</h3>
      <div class="tokens">{ref_tokens}</div>
    </div>
    <div class="panel">
      <h2>Hit</h2>
      <div><strong>SMILES:</strong> {html.escape(match.hit.smiles)}</div>
      <div class="structure">{hit_svg}</div>
      <h3>MolCodon sequence</h3>
      <div class="tokens">{hit_tokens}</div>
    </div>
  </div>
  <div class="panel" style="margin-top:20px;">
    <h2>Component Alignment</h2>
    {components_html}
  </div>
</body>
</html>
'''
    output_path.write_text(html_doc, encoding='utf-8')


def main() -> None:
    parser = argparse.ArgumentParser(description='Generate MolCodon HTML comparison report')
    parser.add_argument('--reference', required=True, help='Reference SMILES')
    parser.add_argument('--hit', required=True, help='Hit SMILES')
    parser.add_argument('-o', '--output', required=True, help='HTML output file')
    args = parser.parse_args()
    match = match_smiles(args.reference, args.hit)
    build_report(match, Path(args.output), title='MolCodon visual match report')
    print(f'Wrote HTML report to {args.output}')


if __name__ == '__main__':
    main()
