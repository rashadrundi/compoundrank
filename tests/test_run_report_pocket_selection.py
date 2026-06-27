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


    def test_report_includes_optional_pose_recovery_benchmark(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)

            summary = {
                "reference_ligand": (
                    "/reference/"
                    "8HUR_7YY_ensitrelvir_crystal.sdf"
                ),
                "poses_sdf": (
                    "/work/docking/"
                    "ensitrelvir/poses.sdf"
                ),
                "rmsd_method": (
                    "symmetry-aware complete "
                    "heavy-atom coordinate RMSD "
                    "without translation or rotation"
                ),
                "bond_order_mapping": (
                    "bond order ignored; atom elements, "
                    "connectivity, ring membership, and "
                    "complete atom coverage required"
                ),
                "rmsd_threshold_angstrom": 2.0,
                "mapped_pose_count": 20,
                "mapping_failure_count": 0,
                "top_cnn_pose": {
                    "pose_index": 1,
                    "receptor_conformer_id": (
                        "snapshot_0001"
                    ),
                    "seed": 20260626,
                    "source_pose_number": 1,
                    "pocket_id": "autobox_01",
                    "cnnscore": 0.823432922,
                    "heavy_atom_rmsd": (
                        1.4892266332345574
                    ),
                },
                "best_sampled_pose": {
                    "pose_index": 7,
                    "receptor_conformer_id": (
                        "submitted_receptor"
                    ),
                    "seed": 20260627,
                    "source_pose_number": 4,
                    "pocket_id": "autobox_01",
                    "cnnscore": 0.701234,
                    "heavy_atom_rmsd": (
                        1.201234
                    ),
                },
                "sampling_pass": True,
                "ranking_pass": True,
                "overall": (
                    "cognate_pose_recovery_"
                    "and_ranking_pass"
                ),
                "evaluated_compound": (
                    "ensitrelvir_prepared"
                ),
                "evaluated_pocket_id": (
                    "autobox_01"
                ),
                "normally_selected_receptor_conformer_id": (
                    "snapshot_0001"
                ),
                "evaluated_receptor_conformer_ids": [
                    "snapshot_0001",
                    "snapshot_0002",
                    "submitted_receptor",
                ],
                "evaluated_receptor_conformer_count": 3,
                "evaluation_stage": (
                    "after normal GNINA scoring, "
                    "PoseBusters filtering, and "
                    "pocket selection"
                ),
                "reference_ligand_used_for_posthoc_evaluation": True,
                "reference_ligand_also_supplied_as_autobox_ligand": True,
                "reference_ligand_used_for_box_definition": True,
                "reference_ligand_used_for_pocket_selection": False,
            }

            (
                output_dir
                / "pose_set_recovery_summary.json"
            ).write_text(
                json.dumps(summary),
                encoding="utf-8",
            )

            (
                output_dir
                / "pose_set_recovery_metrics.csv"
            ).write_text(
                "pose_index,cnnscore,"
                "heavy_atom_rmsd\n"
                "1,0.823432922,"
                "1.4892266332345574\n",
                encoding="utf-8",
            )

            (
                output_dir
                / "pose_set_recovery_report.md"
            ).write_text(
                "# Scored Pose-Recovery Report\n",
                encoding="utf-8",
            )

            report_path = write_run_report(
                output_dir=output_dir
            )

            report = report_path.read_text(
                encoding="utf-8"
            )

            self.assertIn(
                "## Cognate Pose-Recovery Benchmark",
                report,
            )
            self.assertIn(
                "8HUR_7YY_ensitrelvir_crystal.sdf",
                report,
            )
            self.assertIn(
                "Chemically mapped poses | 20",
                report,
            )
            self.assertIn(
                "Mapping failures | 0",
                report,
            )
            self.assertIn(
                "Top CNN score | 0.823433",
                report,
            )
            self.assertIn(
                "Top CNN pose RMSD | 1.489 Å",
                report,
            )
            self.assertIn(
                "Best sampled pose RMSD | 1.201 Å",
                report,
            )
            self.assertIn(
                (
                    "Normally selected receptor conformer: "
                    "snapshot_0001"
                ),
                report,
            )
            self.assertIn(
                (
                    "Evaluated receptor conformers: "
                    "snapshot_0001, snapshot_0002, "
                    "submitted_receptor"
                ),
                report,
            )
            self.assertIn(
                "Evaluated receptor conformer count: 3",
                report,
            )
            self.assertIn(
                (
                    "Top CNN receptor conformer | "
                    "snapshot_0001"
                ),
                report,
            )
            self.assertIn(
                "Top CNN seed | 20260626",
                report,
            )
            self.assertIn(
                "Top CNN source pose | 1",
                report,
            )
            self.assertIn(
                "Top CNN pocket | autobox_01",
                report,
            )
            self.assertIn(
                (
                    "Best sampled receptor conformer | "
                    "submitted_receptor"
                ),
                report,
            )
            self.assertIn(
                "Best sampled seed | 20260627",
                report,
            )
            self.assertIn(
                "Best sampled source pose | 4",
                report,
            )
            self.assertIn(
                "Best sampled pocket | autobox_01",
                report,
            )
            self.assertIn(
                "Sampling pass | yes",
                report,
            )
            self.assertIn(
                "Ranking pass | yes",
                report,
            )
            self.assertIn(
                "Evaluated compound: ensitrelvir_prepared",
                report,
            )
            self.assertIn(
                "Evaluated pocket: autobox_01",
                report,
            )
            self.assertIn(
                (
                    "Reference ligand used for post hoc "
                    "RMSD evaluation: yes"
                ),
                report,
            )
            self.assertIn(
                (
                    "Same reference file supplied as the "
                    "GNINA autobox ligand: yes"
                ),
                report,
            )
            self.assertIn(
                (
                    "Reference ligand used to define the "
                    "docking box: yes"
                ),
                report,
            )
            self.assertIn(
                (
                    "Reference ligand used to choose among "
                    "detected pockets: no"
                ),
                report,
            )
            self.assertIn(
                (
                    "cognate_pose_recovery_"
                    "and_ranking_pass"
                ),
                report,
            )

    def test_report_omits_pose_recovery_when_summary_absent(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)

            report_path = write_run_report(
                output_dir=output_dir
            )

            report = report_path.read_text(
                encoding="utf-8"
            )

            self.assertNotIn(
                "## Cognate Pose-Recovery Benchmark",
                report,
            )


if __name__ == "__main__":
    unittest.main()
