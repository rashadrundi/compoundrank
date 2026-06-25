from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.run_report import (
    write_run_report,
)


class RunReportStructureValidationTests(
    unittest.TestCase
):
    def test_report_includes_geometry_results(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)

            receptor_dir = (
                output
                / "structure_validation"
                / "receptor"
            )

            reference_dir = (
                output
                / "automatic_reference_evidence"
                / "structure_validation"
                / "reference"
            )

            receptor_dir.mkdir(
                parents=True
            )

            reference_dir.mkdir(
                parents=True
            )

            receptor_report = {
                "status": "complete",
                "selection_mode": "report_only",
                "requested_chain": "B",
                "evaluable_residues": 383,
                "summary": {
                    "favored": 346,
                    "allowed": 29,
                    "outliers": 8,
                    "favored_fraction": 346 / 383,
                    "allowed_fraction": 29 / 383,
                    "outlier_fraction": 8 / 383,
                    "screening_flag": (
                        "high_outlier_fraction"
                    ),
                },
            }

            reference_report = {
                "status": "complete",
                "selection_mode": "report_only",
                "requested_chain": "A",
                "evaluable_residues": 380,
                "summary": {
                    "favored": 372,
                    "allowed": 7,
                    "outliers": 1,
                    "favored_fraction": 372 / 380,
                    "allowed_fraction": 7 / 380,
                    "outlier_fraction": 1 / 380,
                    "screening_flag": (
                        "elevated_outlier_fraction"
                    ),
                },
            }

            (
                receptor_dir
                / "ramachandran_validation.json"
            ).write_text(
                json.dumps(
                    receptor_report
                ),
                encoding="utf-8",
            )

            (
                reference_dir
                / "ramachandran_validation.json"
            ).write_text(
                json.dumps(
                    reference_report
                ),
                encoding="utf-8",
            )

            (
                receptor_dir
                / "ramachandran_residues.csv"
            ).write_text(
                "residue,classification\n",
                encoding="utf-8",
            )

            report_path = write_run_report(
                output_dir=output
            )

            report_text = report_path.read_text(
                encoding="utf-8"
            )

            self.assertIn(
                "## Structure Geometry Validation",
                report_text,
            )

            self.assertIn(
                "Submitted receptor",
                report_text,
            )

            self.assertIn(
                "Automatic reference structure",
                report_text,
            )

            self.assertIn(
                "90.34%",
                report_text,
            )

            self.assertIn(
                "high_outlier_fraction",
                report_text,
            )

            self.assertIn(
                "Geometry Review Notes",
                report_text,
            )


if __name__ == "__main__":
    unittest.main()
