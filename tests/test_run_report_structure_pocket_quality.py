from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.run_report import (
    _render_structure_pocket_quality_section,
)


class RunReportStructurePocketQualityTests(
    unittest.TestCase
):
    def test_missing_report_is_omitted(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            lines = (
                _render_structure_pocket_quality_section(
                    Path(temporary)
                )
            )

        self.assertEqual(lines, [])

    def test_renders_global_caution_verdict(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)

            report = {
                "status": "complete",
                "verdict": (
                    "usable_with_global_geometry_caution"
                ),
                "near_box_threshold_angstrom": 4.0,
                "selected_pocket_ids": [
                    "fpocket_02_pocket_2"
                ],
                "outlier_count": 7,
                "inside_selected_box_outliers": [],
                "near_selected_box_outliers": [],
                "selected_box_local_outliers": [],
                "global_ramachandran_summary": {
                    "screening_flag": (
                        "elevated_outlier_fraction"
                    )
                },
            }

            (
                output_dir
                / "structure_pocket_quality.json"
            ).write_text(
                json.dumps(report),
                encoding="utf-8",
            )

            rendered = "\n".join(
                _render_structure_pocket_quality_section(
                    output_dir
                )
            )

        self.assertIn(
            "### Pocket-Localized Structure Quality",
            rendered,
        )
        self.assertIn(
            "Usable with global geometry caution",
            rendered,
        )
        self.assertIn(
            "`fpocket_02_pocket_2`",
            rendered,
        )
        self.assertIn(
            "Global Ramachandran outliers: 7",
            rendered,
        )
        self.assertIn(
            "boxes: 0",
            rendered,
        )
        self.assertIn(
            "near-box threshold: 0",
            rendered,
        )
        self.assertIn(
            "4.0 Å",
            rendered,
        )

    def test_renders_local_geometry_concern(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)

            report = {
                "status": "complete",
                "verdict": (
                    "selected_pocket_geometry_concern"
                ),
                "near_box_threshold_angstrom": 4.0,
                "selected_pocket_ids": [
                    "fpocket_01_pocket_1"
                ],
                "outlier_count": 2,
                "inside_selected_box_outliers": [
                    "ALA:A:10"
                ],
                "near_selected_box_outliers": [],
                "selected_box_local_outliers": [
                    "ALA:A:10"
                ],
                "global_ramachandran_summary": {
                    "screening_flag": (
                        "elevated_outlier_fraction"
                    )
                },
            }

            (
                output_dir
                / "structure_pocket_quality.json"
            ).write_text(
                json.dumps(report),
                encoding="utf-8",
            )

            rendered = "\n".join(
                _render_structure_pocket_quality_section(
                    output_dir
                )
            )

        self.assertIn(
            "Selected pocket geometry concern",
            rendered,
        )
        self.assertIn(
            "inside a selected docking box",
            rendered,
        )
        self.assertIn(
            "`ALA:A:10`",
            rendered,
        )


if __name__ == "__main__":
    unittest.main()
