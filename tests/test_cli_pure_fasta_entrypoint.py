from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.cli import build_parser, main


class CliPureFastaEntrypointTests(unittest.TestCase):
    def test_receptor_is_optional_in_parser(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data_root = root / "data"
            data_root.mkdir()
            fasta = root / "input.faa"
            fasta.write_text(">target\nMSTNPKPQR\n", encoding="utf-8")

            args = build_parser().parse_args(
                [
                    "--fasta",
                    str(fasta),
                    "--data-root",
                    str(data_root),
                    "--output-dir",
                    str(root / "out"),
                    "--skip-structure-prediction",
                ]
            )

            self.assertIsNone(args.receptor)
            self.assertEqual(args.fasta, fasta)

    def test_pure_fasta_without_prediction_writes_clean_skip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data_root = root / "data"
            data_root.mkdir()
            output_dir = root / "out"
            fasta = root / "input.faa"
            fasta.write_text(">target\nMSTNPKPQR\n", encoding="utf-8")

            result = main(
                [
                    "--fasta",
                    str(fasta),
                    "--data-root",
                    str(data_root),
                    "--output-dir",
                    str(output_dir),
                    "--skip-structure-prediction",
                    "--overwrite",
                ]
            )

            self.assertEqual(result, 0)

            skip_path = output_dir / "docking_skipped.json"
            status_path = output_dir / "pure_fasta_structure_status.json"
            report_path = output_dir / "compoundrank_run_report.md"

            self.assertTrue(skip_path.exists())
            self.assertTrue(status_path.exists())
            self.assertTrue(report_path.exists())

            payload = json.loads(skip_path.read_text(encoding="utf-8"))
            self.assertEqual(
                payload["reason_code"],
                "no_receptor_structure_available",
            )
            self.assertEqual(
                payload["pipeline_outcome"],
                "completed_without_docking",
            )
            self.assertEqual(
                payload["downstream_stages"]["gnina_docking"],
                "skipped",
            )


if __name__ == "__main__":
    unittest.main()
