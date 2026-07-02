import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.pipeline import (
    _read_receptor_ramachandran_structure_quality_failure,
)


class PipelineRamachandranHardStopTests(unittest.TestCase):
    def _write_report(
        self,
        root: Path,
        report: dict,
    ) -> Path:
        path = (
            root
            / "structure_validation"
            / "receptor"
            / "ramachandran_validation.json"
        )
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        path.write_text(
            json.dumps(
                report,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return path

    def test_high_outlier_fraction_triggers_structure_quality_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            report_path = self._write_report(
                root,
                {
                    "status": "complete",
                    "selection_mode": "report_only",
                    "source_structure": "/tmp/receptor.pdb",
                    "structure_name": "receptor",
                    "evaluable_residues": 373,
                    "total_polymer_residues": 375,
                    "summary": {
                        "favored_fraction": 0.6621983914209115,
                        "outlier_fraction": 0.19839142091152814,
                        "favored_goal_met": False,
                        "outlier_goal_met": False,
                        "outliers": 74,
                        "screening_flag": "high_outlier_fraction",
                    },
                },
            )

            failure = (
                _read_receptor_ramachandran_structure_quality_failure(
                    root
                )
            )

            self.assertIsNotNone(failure)
            assert failure is not None
            self.assertEqual(
                failure["reason_code"],
                "ramachandran_structure_quality_failed",
            )
            self.assertEqual(
                failure["screening_flag"],
                "high_outlier_fraction",
            )
            self.assertEqual(
                failure["report"],
                str(report_path),
            )
            self.assertAlmostEqual(
                failure["outlier_fraction"],
                0.19839142091152814,
            )

    def test_meets_ramalyze_goals_does_not_trigger_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            self._write_report(
                root,
                {
                    "status": "complete",
                    "selection_mode": "report_only",
                    "summary": {
                        "favored_fraction": 0.991,
                        "outlier_fraction": 0.001,
                        "favored_goal_met": True,
                        "outlier_goal_met": True,
                        "screening_flag": "meets_ramalyze_goals",
                    },
                },
            )

            self.assertIsNone(
                _read_receptor_ramachandran_structure_quality_failure(
                    root
                )
            )

    def test_missing_report_does_not_trigger_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self.assertIsNone(
                _read_receptor_ramachandran_structure_quality_failure(
                    Path(tmpdir)
                )
            )


if __name__ == "__main__":
    unittest.main()
