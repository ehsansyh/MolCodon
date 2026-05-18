#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 -m pip install --upgrade pip
python3 -m pip install "$SCRIPT_DIR"

cat <<'MSG'

MolCodon installed.

Try:
  molcodon-encode "CC(=O)Oc1ccccc1C(=O)O"
  molcodon-decode "SCC CCC CCX NCC NCX CCO CCX OXN SSS"
  molcodon-similarity examples.csv --reference "CCO" -o molcodon_similarity_out
MSG
