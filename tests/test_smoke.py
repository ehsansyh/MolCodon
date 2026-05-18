import subprocess
import sys
import unittest

from rdkit import Chem

from molcodon_decoder import decode_sequence
from molcodon_encoder import encode
from molcodon_match import match_smiles


class MolCodonSmokeTests(unittest.TestCase):
    def test_encode_decode_roundtrip_preserves_inchikey(self):
        smiles = "C[C@H](O)F"
        codons = encode(smiles)
        decoded = decode_sequence(" ".join(codons))
        decoded_smiles = Chem.MolToSmiles(decoded, isomericSmiles=True, canonical=True)

        self.assertEqual(
            Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)),
            Chem.MolToInchiKey(Chem.MolFromSmiles(decoded_smiles)),
        )

    def test_similarity_api_returns_scores(self):
        result = match_smiles("CCO", "CCN")

        self.assertIn("overall", result.scores)
        self.assertGreater(result.scores["overall"], 0)

    def test_encoder_cli_prints_codons(self):
        expected = " ".join(encode("CCO"))
        completed = subprocess.run(
            [sys.executable, "-m", "molcodon_encoder", "CCO"],
            check=True,
            text=True,
            capture_output=True,
        )

        self.assertEqual(completed.stdout.strip(), expected)


if __name__ == "__main__":
    unittest.main()
