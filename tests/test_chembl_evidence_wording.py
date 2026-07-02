from pathlib import Path
import unittest


class ChemblEvidenceWordingTests(unittest.TestCase):
    def test_chembl_reference_wording_does_not_claim_ligand_transfer(self):
        source = Path("compoundrank/chembl_search.py").read_text(
            encoding="utf-8"
        ).lower()

        self.assertNotIn("transferred as ligand evidence", source)
        self.assertNotIn("reference source of transferable ligand evidence", source)
        self.assertIn("manual review", source)
        self.assertIn("does not establish ligand", source)


if __name__ == "__main__":
    unittest.main()
