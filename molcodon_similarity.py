#!/usr/bin/env python3
"""
MolCodon Component Match Scorer
===============================
A completely transparent, explainable scoring system for MolCodon.
Replaces opaque SW/NW alignments and arbitrary weights with exact 
structural component counting (Rings, Branches, Backbone, Bonds).

Usage:
  python molcodon_blast_new.py <csv_file> --reference <SMILES> [--top-n N] [-o outdir]
"""

import argparse
import csv
import html
import time
from collections import Counter
from pathlib import Path

from rdkit import Chem

from molcodon_encoder import encode, EncodeError
from molcodon_decoder import decode_sequence, DecodeError
from molcodon_match import match_smiles
from molcodon_viz import build_report

ATOM_CODON_SET = set(['CCC', 'CCN', 'CCO', 'CCS', 'CNC', 'CNN', 'CNO', 'CNS', 'COC', 'CON'])
BOND_CODON_SET = set(['NCC', 'NCN', 'NCO', 'NCS'])
BRANCH_OPEN_SET = set(['NNC', 'NOC', 'NOS', 'NSC'])
BRANCH_CLOSE_SET = set(['NNN', 'NON', 'NOO', 'NSN'])
RING_OPEN_SET = set(['NNO', 'NSO', 'OSO', 'OCN', 'OCO', 'OCS', 'ONC', 'ONN'])
RING_CLOSE_SET = set(['NNS', 'NSS', 'OSS', 'ONO', 'ONS', 'OSC', 'OSN', 'SCN'])
POSITION_CODON_SET = set([
    'COO', 'COS', 'CSC', 'CSN', 'CSO', 'CSS',
    'OOC', 'OON', 'OOO', 'OOS', 'SCO', 'SCS',
    'SNC', 'SNN', 'SNO', 'SNS',
])
RING_REF_CODON_SET = set(['SOC', 'SON', 'SOO', 'SOS', 'SSC', 'SSN', 'SSO'])
START_CODON = 'SCC'
END_CODON = 'SSS'
FUSION_MARKER = 'OCC'


# ============================================================
# COMPONENT EXTRACTION
# ============================================================

def annotate_codons(codons):
    rows = []
    ring_stack = []
    branch_stack = []
    ring_id_counter = 0
    branch_id_counter = 0
    ring_open_to_close = {
        'NNO': 'NNS', 'NSO': 'NSS', 'OSO': 'OSS', 'OCN': 'ONO',
        'OCO': 'ONS', 'OCS': 'OSC', 'ONC': 'OSN', 'ONN': 'SCN',
    }
    branch_open_to_close = {
        'NNC': 'NNN', 'NOC': 'NON', 'NOS': 'NOO', 'NSC': 'NSN',
    }
    for codon in codons:
        row = {'codon': codon, 'class': 'other', 'ring_id': None, 'branch_id': None}
        if codon in ATOM_CODON_SET:
            row['class'] = 'atom'
        elif codon in BOND_CODON_SET:
            row['class'] = 'bond'
        elif codon == START_CODON:
            row['class'] = 'start'
        elif codon == END_CODON:
            row['class'] = 'end'
        elif codon == FUSION_MARKER:
            row['class'] = 'fusion'
        elif codon in RING_OPEN_SET:
            row['class'] = 'ring:open'
            row['ring_id'] = ring_id_counter
            ring_stack.append((codon, ring_id_counter))
            ring_id_counter += 1
        elif codon in RING_CLOSE_SET:
            row['class'] = 'ring:close'
            if ring_stack:
                _, rid = ring_stack.pop()
                row['ring_id'] = rid
        elif codon in BRANCH_OPEN_SET:
            row['class'] = 'branch:open'
            row['branch_id'] = branch_id_counter
            branch_stack.append((codon, branch_id_counter))
            branch_id_counter += 1
        elif codon in BRANCH_CLOSE_SET:
            row['class'] = 'branch:close'
            if branch_stack:
                _, bid = branch_stack.pop()
                row['branch_id'] = bid
        rows.append(row)
    return rows


def extract_components(codons, rows):
    """
    Extract exact components from a MolCodon sequence:
    1. rings: list of tuple(codons inside ring)
    2. branches: list of tuple(codons inside branch)
    3. backbone: list of (atom and bond codons not in ring/branch)
    4. bond_types: multiset of all bond codons
    """
    rings = []
    branches = []
    backbone = []
    bond_types = Counter()
    
    # We use the annotated rows to cleanly group tokens
    ring_stack = []
    branch_stack = []
    
    for row in rows:
        c = row['codon']
        cls = row['class']
        
        # Track bond types globally for composition
        if c in BOND_CODON_SET:
            bond_types[c] += 1
            
        # Ring extraction
        if cls == 'ring:open':
            ring_stack.append({'id': row['ring_id'], 'codons': []})
            continue
        elif cls == 'ring:close':
            if ring_stack:
                r = ring_stack.pop()
                # Store as tuple so it's hashable/comparable
                rings.append(tuple(r['codons']))
            continue
            
        # Branch extraction
        if cls == 'branch:open':
            branch_stack.append({'id': row['branch_id'], 'codons': []})
            continue
        elif cls == 'branch:close':
            if branch_stack:
                b = branch_stack.pop()
                branches.append(tuple(b['codons']))
            continue
            
        # If we are inside a ring or branch, add to the innermost one
        if branch_stack:
            branch_stack[-1]['codons'].append(c)
        elif ring_stack:
            ring_stack[-1]['codons'].append(c)
        else:
            # Main backbone (exclude control/position markers if desired, 
            # but keeping them ensures strict topological matching)
            if c not in (START_CODON, END_CODON):
                backbone.append(c)
                
    return {
        'rings': Counter(rings),
        'branches': Counter(branches),
        'backbone': tuple(backbone),
        'bond_types': bond_types
    }

# ============================================================
# MAIN PIPELINE
# ============================================================

def analyze(ref_smiles: str, candidates: list):
    ref_codons = encode(ref_smiles)
    ref_mol = decode_sequence(' '.join(ref_codons))
    ref_canonical = Chem.MolToSmiles(ref_mol, canonical=True)

    ref = {
        'smiles': ref_smiles,
        'canonical': ref_canonical,
        'codons': ref_codons,
        'meta': {'canonical': ref_canonical},
    }
    
    results = []
    for cand in candidates:
        smi = cand['smiles']
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            decoded = decode_sequence(' '.join(encode(smi)))
            canonical = Chem.MolToSmiles(decoded, canonical=True)
            match = match_smiles(ref_canonical, canonical)
            
            results.append({
                'name': cand['name'],
                'smiles': smi,
                'canonical': canonical,
                'codons': match.hit.tokens,
                'meta': {'canonical': canonical},
                'scores': match.scores,
                'tanimoto': match.tanimoto,
                'match': match,
            })
        except (EncodeError, DecodeError, Exception):
            continue
            
    results.sort(key=lambda r: (r['scores']['overall'], r['tanimoto']), reverse=True)
    return ref, results

def write_results_csv(path: Path, results: list):
    fields = [
        'rank', 'name', 'smiles', 'canonical', 'num_codons',
        'overall_score', 'ring_match', 'ring_f1', 'branch_match', 'branch_f1',
        'attachment_match', 'attachment_f1',
        'bond_match', 'bond_f1', 'pharmacophore_match', 'pharmacophore_f1', 'backbone_score', 'tanimoto'
    ]
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rank, row in enumerate(results, start=1):
            sc = row['scores']
            w.writerow({
                'rank': rank, 'name': row['name'],
                'smiles': row['smiles'], 'canonical': row['canonical'],
                'num_codons': len(row['codons']),
                'overall_score': sc['overall'],
                'ring_match': sc['ring_match'],
                'ring_f1': sc['ring_f1'],
                'branch_match': sc['branch_match'],
                'branch_f1': sc['branch_f1'],
                'attachment_match': sc.get('attachment_match', ''),
                'attachment_f1': sc.get('attachment_f1', ''),
                'bond_match': sc['bond_match'],
                'bond_f1': sc['bond_f1'],
                'pharmacophore_match': sc.get('pharmacophore_match', ''),
                'pharmacophore_f1': sc.get('pharmacophore_f1', ''),
                'backbone_score': sc['backbone'],
                'tanimoto': row['tanimoto']
            })

def _build_index_table(results_subset: list, link_prefix: str) -> str:
    rows = []
    for idx, r in enumerate(results_subset, start=1):
        sc = r['scores']
        link = f'{link_prefix}_{idx:02d}.html'
        rows.append(
            f'<tr><td>{idx}</td>'
            f'<td><a href="{link}">{html.escape(r["name"])}</a></td>'
            f'<td>{sc["overall"]:.2f}</td>'
            f'<td>{sc["ring_f1"]:.1f} ({sc["ring_match"]})</td>'
            f'<td>{sc["branch_f1"]:.1f} ({sc["branch_match"]})</td>'
            f'<td>{sc["attachment_f1"]:.1f} ({sc["attachment_match"]})</td>'
            f'<td>{sc["pharmacophore_f1"]:.1f}</td>'
            f'<td>{r["tanimoto"]:.2f}</td></tr>'
        )
    return ''.join(rows)


def build_index_html(path: Path, ref_canonical: str, results: list,
                     molcodon_count: int, tanimoto_count: int):
    INDEX_STYLE = '''<style>
body { font-family: Inter, Arial, sans-serif; margin: 24px; color: #0f172a; background: #f8fafc; }
.panel { background: white; border-radius: 16px; padding: 24px; box-shadow: 0 8px 24px rgba(15,23,42,0.08); max-width: 1100px; margin: 0 auto 24px; }
h1,h2 { margin-top: 0; }
table { width: 100%; border-collapse: collapse; margin-top: 12px; }
th, td { border-bottom: 1px solid #e2e8f0; padding: 10px 12px; text-align: left; font-size: 14px; }
th { color: #334155; font-weight: 600; background: #f8fafc; }
a { color: #2563eb; text-decoration: none; font-weight: 500; }
a:hover { text-decoration: underline; }
.reference { background: #f1f5f9; padding: 16px; border-radius: 8px; margin-bottom: 20px; font-family: monospace; word-break: break-all; font-size: 13px; }
</style>'''

    HEADER = '<tr><th>#</th><th>Name</th><th>Overall</th><th>Ring</th><th>Branch</th><th>Attachment</th><th>Pharm</th><th>Tanimoto</th></tr>'

    # MolCodon-ranked table (results already sorted by overall)
    molcodon_rows = _build_index_table(results[:molcodon_count], 'report_rank')

    # Tanimoto-ranked table
    tani_sorted = sorted(results, key=lambda r: r['tanimoto'], reverse=True)
    tani_rows = _build_index_table(tani_sorted[:tanimoto_count], 'report_tani')

    html_content = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"/><title>MolCodon BLAST Results</title>{INDEX_STYLE}</head>
<body>
<div class="panel">
  <h1>MolCodon BLAST Results</h1>
  <div class="reference"><strong>Reference:</strong> {html.escape(ref_canonical)}</div>
</div>
<div class="panel">
  <h2>Top by MolCodon Overall Score</h2>
  <table><thead>{HEADER}</thead><tbody>{molcodon_rows}</tbody></table>
</div>
<div class="panel">
  <h2>Top by Tanimoto Similarity</h2>
  <table><thead>{HEADER}</thead><tbody>{tani_rows}</tbody></table>
</div>
</body></html>'''
    path.write_text(html_content, encoding='utf-8')

def main():
    t0 = time.time()
    parser = argparse.ArgumentParser(description='MolCodon Component Match')
    parser.add_argument('csvfile', help='CSV file with SMILES column')
    parser.add_argument('--reference', required=True, help='Reference SMILES string')
    parser.add_argument('--output-dir', '-o', default='molcodon_blast_out', help='Output directory')
    parser.add_argument('--html-top', type=int, default=1, help='Generate HTML reports for top N hits')
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Simple CSV loader
    candidates = []
    with open(args.csvfile) as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, 1):
            # Guess SMILES and Name columns
            smi_col = next((c for c in row if c.lower() in ('smiles','smi')), None)
            name_col = next((c for c in row if c.lower() in ('name','id','compound')), None)
            
            smi = row.get(smi_col, '').strip() if smi_col else list(row.values())[0].strip()
            if not smi: continue
            name = row.get(name_col, '').strip() if name_col else f'mol_{idx}'
            candidates.append({'name': name, 'smiles': smi})
            
    print(f'Loaded {len(candidates)} candidates.')
    
    ref, results = analyze(args.reference, candidates)
    
    write_results_csv(outdir / 'results_component.csv', results)

    # ── Generate HTML reports for top MolCodon-ranked hits ───
    molcodon_written = 0
    for idx, result in enumerate(results[:max(0, args.html_top)], start=1):
        try:
            match = result.get('match') or match_smiles(ref['canonical'], result['canonical'])
            build_report(match, outdir / f'report_rank_{idx:02d}.html',
                         title=f'MolCodon #{idx}: {result["name"]}')
            molcodon_written += 1
        except Exception as e:
            print(f"Warning: report failed for {result['name']}: {e}")

    # ── Generate HTML reports for top Tanimoto-ranked hits ──
    tani_sorted = sorted(results, key=lambda r: r['tanimoto'], reverse=True)
    tani_written = 0
    for idx, result in enumerate(tani_sorted[:max(0, args.html_top)], start=1):
        try:
            match = result.get('match') or match_smiles(ref['canonical'], result['canonical'])
            build_report(match, outdir / f'report_tani_{idx:02d}.html',
                         title=f'Tanimoto #{idx}: {result["name"]}')
            tani_written += 1
        except Exception as e:
            print(f"Warning: tanimoto report failed for {result['name']}: {e}")

    # ── Build index page ───────────────────────────────────
    if molcodon_written > 0 or tani_written > 0:
        build_index_html(outdir / 'index.html', ref['canonical'], results,
                         molcodon_written, tani_written)

    runtime = time.time() - t0
    print(f'\nTop 10 by MolCodon Overall for {ref["canonical"]}:')
    print(f'{"Rank":<5} | {"Name":<20} | {"Overall":<8} | {"Rings":<8} | {"Branch":<8} | {"Attach":<8} | {"Tanimoto":<8}')
    print('-' * 80)
    for i, r in enumerate(results[:10], 1):
        sc = r['scores']
        print(f'{i:<5} | {r["name"][:20]:<20} | {sc["overall"]:<8.2f} | {sc["ring_match"]:<8} | {sc["branch_match"]:<8} | {sc["attachment_match"]:<8} | {r["tanimoto"]:<8.2f}')

    print(f'\nTop 10 by Tanimoto:')
    for i, r in enumerate(tani_sorted[:10], 1):
        sc = r['scores']
        print(f'{i:<5} | {r["name"][:20]:<20} | {sc["overall"]:<8.2f} | {r["tanimoto"]:<8.2f}')

    print(f'\nResults: {outdir / "results_component.csv"}')
    print(f'HTML reports: {molcodon_written} MolCodon + {tani_written} Tanimoto')
    print(f'Runtime: {runtime:.2f}s')

if __name__ == '__main__':
    main()
