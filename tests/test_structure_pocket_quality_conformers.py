from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.structure_pocket_quality import (
    evaluate_structure_pocket_quality,
    run_structure_pocket_quality,
)


class StructurePocketQualityConformerTests(
    unittest.TestCase
):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

        self.receptor = self.root / "receptor.pdb"
        self.receptor.write_text(
            (
                "ATOM      1  N   ALA A   1       "
                "0.000   0.000   0.000  1.00 20.00"
                "           N  \n"
                "END\n"
            ),
            encoding="utf-8",
        )

        self.ramachandran = self._write_json(
            "ramachandran.json",
            {
                "status": "complete",
                "model_index": 0,
                "evaluable_residues": 1,
                "summary": {
                    "favored": 1,
                    "favored_fraction": 1.0,
                    "allowed": 0,
                    "allowed_fraction": 0.0,
                    "outliers": 0,
                    "outlier_fraction": 0.0,
                    "screening_flag": (
                        "meets_ramalyze_goals"
                    ),
                },
                "goals": {
                    "favored_fraction_goal_met": True,
                    "outlier_fraction_goal_met": True,
                },
                "residues": [],
            },
        )

        self.pockets = self._write_json(
            "pocket_definitions.json",
            {
                "pockets": [
                    {
                        "pocket_id": "pocket_1",
                        "center_x": 0.0,
                        "center_y": 0.0,
                        "center_z": 0.0,
                        "size_x": 10.0,
                        "size_y": 10.0,
                        "size_z": 10.0,
                    },
                    {
                        "pocket_id": "pocket_2",
                        "center_x": 20.0,
                        "center_y": 20.0,
                        "center_z": 20.0,
                        "size_x": 10.0,
                        "size_y": 10.0,
                        "size_z": 10.0,
                    },
                ]
            },
        )

        self.selection = self._write_json(
            "pocket_selection_summary.json",
            {
                "selected_pockets": [
                    {
                        "selected": True,
                        "compound": "compound_1",
                        "pocket_id": "pocket_1",
                        "receptor_conformer_id": (
                            "snapshot_0001"
                        ),
                    },
                    {
                        "selected": True,
                        "compound": "compound_2",
                        "pocket_id": "pocket_2",
                        "receptor_conformer_id": (
                            "snapshot_0002"
                        ),
                    },
                ]
            },
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write_json(
        self,
        filename: str,
        payload: dict,
    ) -> Path:
        path = self.root / filename
        path.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )
        return path

    def test_filters_selected_rows_by_conformer(
        self,
    ) -> None:
        report = evaluate_structure_pocket_quality(
            structure_path=self.receptor,
            ramachandran_report_path=(
                self.ramachandran
            ),
            pocket_definitions_path=self.pockets,
            pocket_selection_summary_path=(
                self.selection
            ),
            receptor_conformer_id=(
                "snapshot_0002"
            ),
        )

        self.assertEqual(
            report["receptor_conformer_id"],
            "snapshot_0002",
        )
        self.assertEqual(
            report["selected_pocket_ids"],
            ["pocket_2"],
        )

    def test_missing_selected_conformer_raises(
        self,
    ) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "snapshot_0003",
        ):
            evaluate_structure_pocket_quality(
                structure_path=self.receptor,
                ramachandran_report_path=(
                    self.ramachandran
                ),
                pocket_definitions_path=(
                    self.pockets
                ),
                pocket_selection_summary_path=(
                    self.selection
                ),
                receptor_conformer_id=(
                    "snapshot_0003"
                ),
            )

    def test_runner_uses_conformer_directory(
        self,
    ) -> None:
        output_dir = self.root / "output"

        result = run_structure_pocket_quality(
            structure_path=self.receptor,
            ramachandran_report_path=(
                self.ramachandran
            ),
            pocket_definitions_path=self.pockets,
            pocket_selection_summary_path=(
                self.selection
            ),
            output_dir=output_dir,
            receptor_conformer_id=(
                "snapshot_0002"
            ),
        )

        expected = (
            output_dir
            / "structure_pocket_quality"
            / "snapshot_0002"
            / "structure_pocket_quality.json"
        )

        self.assertEqual(
            result["output_path"],
            expected,
        )
        self.assertTrue(expected.is_file())

        saved = json.loads(
            expected.read_text(encoding="utf-8")
        )

        self.assertEqual(
            saved["receptor_conformer_id"],
            "snapshot_0002",
        )
        self.assertEqual(
            saved["selected_pocket_ids"],
            ["pocket_2"],
        )


if __name__ == "__main__":
    unittest.main()
