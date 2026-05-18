# MolCodon

MolCodon is a molecular representation toolkit for encoding SMILES strings into codon-like chemical token sequences, decoding MolCodon sequences back into RDKit molecules, and scoring component-level molecular similarity.

The package includes:

- `molcodon_encoder.py`: SMILES to MolCodon sequence.
- `molcodon_decoder.py`: MolCodon sequence to RDKit molecule or canonical SMILES.
- `molcodon_match.py`: explainable component-based pairwise similarity.
- `molcodon_similarity.py`: CSV search workflow with MolCodon and Tanimoto rankings.
- `molcodon_viz.py`: optional HTML reports for matched molecular components.
- `MOLCODON_CODON_DICTIONARY.md`: codon dictionary.

## Requirements

MolCodon requires Python 3.9+ and RDKit.

Recommended RDKit installation:

```bash
conda install -c conda-forge rdkit
```

Pip installation can also work in many environments:

```bash
python3 -m pip install rdkit
```

## Install

From the repository directory:

```bash
chmod +x install.sh
./install.sh
```

Or install directly with pip:

```bash
python3 -m pip install .
```

For development:

```bash
python3 -m pip install -e .
```

## Quick Start

Encode a SMILES string:

```bash
molcodon-encode "CC(=O)Oc1ccccc1C(=O)O"
```

Decode a MolCodon sequence:

```bash
molcodon-decode "SCC CCC CCX NCC NCX CCO CCX OXN SSS"
```

Use the Python API:

```python
from rdkit import Chem
from molcodon_encoder import encode
from molcodon_decoder import decode_sequence

smiles = "C[C@H](O)F"
codons = encode(smiles)
decoded = decode_sequence(" ".join(codons))
decoded_smiles = Chem.MolToSmiles(decoded, isomericSmiles=True, canonical=True)

print(codons)
print(decoded_smiles)
```

## Similarity Search

Prepare a CSV file with a `smiles` or `smi` column. A `name`, `id`, or `compound` column is optional.

```csv
name,smiles
ethanol,CCO
ethylamine,CCN
aspirin,CC(=O)Oc1ccccc1C(=O)O
```

Run MolCodon similarity search:

```bash
molcodon-similarity examples.csv --reference "CCO" -o molcodon_similarity_out --html-top 3
```

Outputs:

- `results_component.csv`: ranked similarity table.
- `index.html`: summary page when HTML reports are requested.
- `report_rank_*.html`: top MolCodon-ranked component reports.
- `report_tani_*.html`: top Tanimoto-ranked reports.

## Similarity Scores

MolCodon similarity is component based:

- `overall`: weighted component-level MolCodon similarity.
- `ring_f1`: overlap of normalized ring components.
- `branch_f1`: overlap of branch components.
- `attachment_f1`: agreement of branch attachment context.
- `bond_f1`: bond type overlap.
- `pharmacophore_f1`: overlap of pharmacophore annotations.
- `backbone`: token n-gram similarity of non-ring/non-branch backbone.
- `tanimoto`: Morgan fingerprint Tanimoto similarity from RDKit.

## Notes

MolCodon currently supports common organic atoms encoded by the included codon dictionary: C, N, O, S, F, Cl, Br, I, P, and B. Unsupported atoms or unsupported graph patterns raise an `EncodeError` or `DecodeError`.

The encoder and decoder preserve standard atom chirality (`R/S`) and double-bond stereochemistry (`E/Z`) for supported molecules. Pseudoasymmetric lowercase CIP descriptors (`r/s`) are a known limitation unless explicitly extended with additional codons.

## Test

Run the smoke tests from the repository root:

```bash
PYTHONPATH=. python3 -m unittest discover -s tests
```

## Citation

If you use MolCodon in a manuscript, cite the associated MolCodon paper or preprint when available.
