import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from compoundrank.cli import build_parser, main


class CliAutoRetrieveTests(unittest.TestCase):
    def test_parser_accepts_auto_retrieve_flags(self):
        args = build_parser().parse_args(
            [
                "--receptor", "/tmp/receptor.pdb",
                "--data-root", "/tmp/data",
                "--fasta", "/tmp/protein.faa",
                "--auto-retrieve-ligands",
                "--auto-retrieve-max-candidates", "4",
                "--auto-retrieve-no-fetch-structures",
            ]
        )

        self.assertTrue(args.auto_retrieve_ligands)
        self.assertEqual(args.auto_retrieve_max_candidates, 4)
        self.assertTrue(args.auto_retrieve_no_fetch_structures)

    def test_auto_retrieve_requires_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            receptor = tmp / "receptor.pdb"
            receptor.write_text("END\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                main(
                    [
                        "--receptor", str(receptor),
                        "--data-root", str(tmp / "data"),
                        "--auto-retrieve-ligands",
                    ]
                )

    def test_auto_retrieve_passes_empty_manual_requests_to_pipeline(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            receptor = tmp / "receptor.pdb"
            fasta = tmp / "protein.faa"
            receptor.write_text("END\n", encoding="utf-8")
            fasta.write_text(">protein\nAAAA\n", encoding="utf-8")

            with patch("compoundrank.cli.run_pipeline") as mocked:
                main(
                    [
                        "--receptor", str(receptor),
                        "--data-root", str(tmp / "data"),
                        "--fasta", str(fasta),
                        "--auto-retrieve-ligands",
                        "--auto-retrieve-max-candidates", "3",
                    ]
                )

            kwargs = mocked.call_args.kwargs
            self.assertEqual(kwargs["ligand_requests"], [])
            self.assertTrue(kwargs["auto_retrieve_ligands"])
            self.assertEqual(kwargs["auto_retrieve_max_candidates"], 3)


if __name__ == "__main__":
    unittest.main()
