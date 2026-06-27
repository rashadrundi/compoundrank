from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.run_report import (
    _render_structure_pocket_quality_section,
)


class RunReportStructurePocketQualityEnsembleTests(
    unittest.TestCase
):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write_json(
        self,
        filename: str,
        payload: dict,
    ) -> Path:
        path = self.root / filename
        path.write_text(
            json.dumps(
                payload,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return path

    def _render(
        self,
    ) -> str:
        return "\n".join(
            _render_structure_pocket_quality_section(
                self.root
            )
        )

    def test_renders_selected_conformer_table(
        self,
    ) -> None:
        conformer_report = (
            self.root
            / "structure_pocket_quality"
            / "snapshot_0002"
            / "structure_pocket_quality.json"
        )
        conformer_report.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        conformer_report.write_text(
            "{}\n",
            encoding="utf-8",
        )

        self._write_json(
            "structure_pocket_quality_ensemble.json",
            {
                "schema_version": (
                    "structure_pocket_quality_"
                    "ensemble.v0.2"
                ),
                "status": "complete",
                "selected_conformer_count": 1,
                "completed_conformer_count": 1,
                "incomplete_conformer_count": 0,
                "overall_verdict": (
                    "usable_with_global_"
                    "geometry_caution"
                ),
                "conformers": [
                    {
                        "receptor_conformer_id": (
                            "snapshot_0002"
                        ),
                        "status": "complete",
                        "selected_pocket_ids": [
                            "autobox_01"
                        ],
                        "outlier_count": 5,
                        "selected_box_local_outliers": [
                            "HIS:B:296"
                        ],
                        "selected_box_local_outlier_count": 1,
                        "selected_pose_local_outliers": [],
                        "selected_pose_local_outlier_count": 0,
                        "box_edge_only_outliers": [
                            "HIS:B:296"
                        ],
                        "box_edge_only_outlier_count": 1,
                        "selected_pose_available": True,
                        "verdict": (
                            "usable_with_global_"
                            "geometry_caution"
                        ),
                        "report_path": str(
                            conformer_report
                        ),
                    }
                ],
            },
        )

        rendered = self._render()

        self.assertIn(
            "Overall verdict: "
            "**Usable with global geometry caution**",
            rendered,
        )
        self.assertIn(
            "`snapshot_0002`",
            rendered,
        )
        self.assertIn(
            "`autobox_01`",
            rendered,
        )
        self.assertIn(
            (
                "| Global outliers | Box-local | "
                "Pose-local | Box-edge only |"
            ),
            rendered,
        )
        self.assertIn(
            "| 5 | 1 | 0 | 1 |",
            rendered,
        )
        self.assertIn(
            (
                "no identified backbone outlier "
                "triggered a selected-pose-local "
                "geometry concern"
            ),
            rendered,
        )
        self.assertIn(
            (
                "Box-local-only advisories may still "
                "be present"
            ),
            rendered,
        )
        self.assertNotIn(
            (
                "no identified backbone outlier was "
                "localized inside or near a selected "
                "docking box"
            ),
            rendered,
        )
        self.assertIn(
            (
                "structure_pocket_quality/"
                "snapshot_0002/"
                "structure_pocket_quality.json"
            ),
            rendered,
        )

    def test_renders_incomplete_conformer_reason(
        self,
    ) -> None:
        self._write_json(
            "structure_pocket_quality_ensemble.json",
            {
                "status": "partial",
                "selected_conformer_count": 2,
                "completed_conformer_count": 1,
                "incomplete_conformer_count": 1,
                "overall_verdict": (
                    "manual_review_incomplete_"
                    "conformer_quality"
                ),
                "conformers": [
                    {
                        "receptor_conformer_id": (
                            "snapshot_0001"
                        ),
                        "status": "skipped",
                        "reason": (
                            "Conformer-specific "
                            "Ramachandran validation "
                            "was not complete."
                        ),
                    }
                ],
            },
        )

        rendered = self._render()

        self.assertIn(
            "**Partial**",
            rendered,
        )
        self.assertIn(
            "`snapshot_0001`",
            rendered,
        )
        self.assertIn(
            "Not evaluated",
            rendered,
        )
        self.assertIn(
            (
                "Conformer-specific Ramachandran "
                "validation was not complete."
            ),
            rendered,
        )

    def test_aggregate_precedes_legacy_report(
        self,
    ) -> None:
        self._write_json(
            "structure_pocket_quality.json",
            {
                "status": "complete",
                "verdict": (
                    "selected_pocket_geometry_concern"
                ),
                "selected_pocket_ids": [
                    "legacy_pocket"
                ],
            },
        )

        self._write_json(
            "structure_pocket_quality_ensemble.json",
            {
                "status": "complete",
                "selected_conformer_count": 1,
                "completed_conformer_count": 1,
                "incomplete_conformer_count": 0,
                "overall_verdict": "strong",
                "conformers": [],
            },
        )

        rendered = self._render()

        self.assertIn(
            "Overall verdict: **Strong**",
            rendered,
        )
        self.assertNotIn(
            "legacy_pocket",
            rendered,
        )


if __name__ == "__main__":
    unittest.main()
