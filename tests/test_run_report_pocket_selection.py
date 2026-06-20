from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.run_report import write_run_report


class RunReportPocketSelectionTests(
    unittest.TestCase
):
    def test_report_includes_selected_and_alternative_pockets(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)

            selection = {
                "selection_method": (
                    "highest top-pose CNNscore"
                ),
                "reference_ligand_used_for_selection": False,
                "compound_count": 1,
                "attempt_count": 3,
                "attempts": [
                    {
                        "compound": "darunavir",
                        "selection_rank": 1,
                        "selected": True,
                        "pocket_id": "fpocket_03_pocket_3",
                        "pocket_rank": 3,
                        "fpocket_score": 0.126,
                        "raw_poses": 20,
                        "accepted_poses": 20,
                        "rejected_poses": 0,
                        "score_source": "accepted_poses",
                        "top_cnn_score": 0.9818345904,
                        "top_cnn_affinity": 10.3630628586,
                        "top_minimized_affinity": -9.42787,
                    },
                    {
                        "compound": "darunavir",
                        "selection_rank": 2,
                        "selected": False,
                        "pocket_id": "fpocket_02_pocket_2",
                        "pocket_rank": 2,
                        "fpocket_score": 0.145,
                        "raw_poses": 20,
                        "accepted_poses": 20,
                        "rejected_poses": 0,
                        "score_source": "accepted_poses",
                        "top_cnn_score": 0.6176877618,
                        "top_cnn_affinity": 7.1356339455,
                        "top_minimized_affinity": -2.0432,
                    },
                    {
                        "compound": "darunavir",
                        "selection_rank": 3,
                        "selected": False,
                        "pocket_id": "fpocket_01_pocket_1",
                        "pocket_rank": 1,
                        "fpocket_score": 0.147,
                        "raw_poses": 20,
                        "accepted_poses": 20,
                        "rejected_poses": 0,
                        "score_source": "accepted_poses",
                        "top_cnn_score": 0.615578413,
                        "top_cnn_affinity": 4.9278774261,
                        "top_minimized_affinity": 0.60846,
                    },
                ],
            }

            definitions = {
                "pocket_count": 3,
                "ranking_method": (
                    "fpocket score descending"
                ),
                "reference_ligand_used_for_selection": False,
                "pockets": [
                    {
                        "pocket_id": "fpocket_03_pocket_3",
                        "pocket_rank": 3,
                        "fpocket_score": 0.126,
                        "center_x": 7.7325,
                        "center_y": -15.8355,
                        "center_z": -2.696,
                        "size_x": 21.731,
                        "size_y": 20.0,
                        "size_z": 24.512,
                        "source": "fpocket rank 3",
                    }
                ],
            }

            (
                output_dir
                / "pocket_selection_summary.json"
            ).write_text(
                json.dumps(selection),
                encoding="utf-8",
            )

            (
                output_dir
                / "pocket_definitions.json"
            ).write_text(
                json.dumps(definitions),
                encoding="utf-8",
            )

            report_path = write_run_report(
                output_dir=output_dir
            )
            report = report_path.read_text(
                encoding="utf-8"
            )

            self.assertIn(
                "## Pocket Detection and Selection",
                report,
            )
            self.assertIn(
                "fpocket_03_pocket_3",
                report,
            )
            self.assertIn(
                "fpocket_02_pocket_2",
                report,
            )
            self.assertIn(
                "0.9818",
                report,
            )
            self.assertIn(
                "0.3641",
                report,
            )
            self.assertIn(
                "high separation",
                report,
            )
            self.assertIn(
                "GNINA selected fpocket rank 3",
                report,
            )
            self.assertIn(
                "Reference ligand used for pocket selection: no",
                report,
            )
            self.assertIn(
                "7.732, -15.835, -2.696",
                report,
            )

    def test_report_warns_about_low_margin_and_raw_fallback(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)

            selection = {
                "selection_method": (
                    "highest top-pose CNNscore"
                ),
                "reference_ligand_used_for_selection": False,
                "compound_count": 1,
                "attempt_count": 2,
                "attempts": [
                    {
                        "compound": "example",
                        "selection_rank": 1,
                        "selected": True,
                        "pocket_id": "pocket_2",
                        "pocket_rank": 2,
                        "raw_poses": 20,
                        "accepted_poses": 0,
                        "rejected_poses": 20,
                        "score_source": "raw_pose_fallback",
                        "top_cnn_score": 0.61,
                        "top_cnn_affinity": 7.0,
                        "top_minimized_affinity": -6.0,
                    },
                    {
                        "compound": "example",
                        "selection_rank": 2,
                        "selected": False,
                        "pocket_id": "pocket_1",
                        "pocket_rank": 1,
                        "raw_poses": 20,
                        "accepted_poses": 20,
                        "rejected_poses": 0,
                        "score_source": "accepted_poses",
                        "top_cnn_score": 0.59,
                        "top_cnn_affinity": 6.9,
                        "top_minimized_affinity": -5.9,
                    },
                ],
            }

            (
                output_dir
                / "pocket_selection_summary.json"
            ).write_text(
                json.dumps(selection),
                encoding="utf-8",
            )

            report_path = write_run_report(
                output_dir=output_dir
            )
            report = report_path.read_text(
                encoding="utf-8"
            )

            self.assertIn(
                "low separation",
                report,
            )
            self.assertIn(
                "raw GNINA poses",
                report,
            )
            self.assertIn(
                "treated as uncertain",
                report,
            )


if __name__ == "__main__":
    unittest.main()
