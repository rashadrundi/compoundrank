from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D

from compoundrank.structure_pocket_quality import (
    evaluate_structure_pocket_quality,
    run_structure_pocket_quality,
)


class StructurePocketQualityTests(
    unittest.TestCase
):
    def setUp(self) -> None:
        self.temporary_directory = (
            tempfile.TemporaryDirectory()
        )
        self.root = Path(
            self.temporary_directory.name
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _write_json(
        self,
        name: str,
        data: dict,
    ) -> Path:
        path = self.root / name
        path.write_text(
            json.dumps(
                data,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return path

    def _write_receptor(
        self,
        *,
        x: float,
        y: float = 0.0,
        z: float = 0.0,
    ) -> Path:
        path = self.root / "receptor.pdb"

        pdb_text = (
            "ATOM      1  N   ALA A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            "  1.00 20.00           N  \n"
            "ATOM      2  CA  ALA A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            "  1.00 20.00           C  \n"
            "ATOM      3  C   ALA A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            "  1.00 20.00           C  \n"
            "TER\n"
            "END\n"
        )

        path.write_text(
            pdb_text,
            encoding="utf-8",
        )

        return path

    def _write_ramachandran_report(
        self,
        *,
        include_outlier: bool = True,
        goals_met: bool = False,
    ) -> Path:
        residues = []

        if include_outlier:
            residues.append(
                {
                    "model": "1",
                    "chain": "A",
                    "residue_name": "ALA",
                    "residue_number": "1",
                    "residue": "ALA:A:1",
                    "category": "general",
                    "classification": "outlier",
                    "score": 0.0001,
                }
            )

        report = {
            "schema_version": (
                "ramachandran_validation.v0.2"
            ),
            "status": "complete",
            "model_index": 0,
            "summary": {
                "favored_goal_met": goals_met,
                "outlier_goal_met": goals_met,
                "favored_fraction": (
                    1.0 if goals_met else 0.95
                ),
                "outlier_fraction": (
                    0.0
                    if goals_met
                    else 0.01
                ),
            },
            "residues": residues,
        }

        return self._write_json(
            "ramachandran_validation.json",
            report,
        )

    def _write_pocket_definitions(
        self,
    ) -> Path:
        return self._write_json(
            "pocket_definitions.json",
            {
                "pocket_count": 1,
                "pockets": [
                    {
                        "pocket_id": (
                            "fpocket_01_pocket_1"
                        ),
                        "pocket_rank": 1,
                        "fpocket_score": 0.5,
                        "mode": "explicit",
                        "center_x": 0.0,
                        "center_y": 0.0,
                        "center_z": 0.0,
                        "size_x": 20.0,
                        "size_y": 20.0,
                        "size_z": 20.0,
                    }
                ],
            },
        )

    def _write_selection_summary(
        self,
    ) -> Path:
        return self._write_json(
            "pocket_selection_summary.json",
            {
                "selected_pockets": [
                    {
                        "selected": True,
                        "compound": "test_compound",
                        "pocket_id": (
                            "fpocket_01_pocket_1"
                        ),
                    }
                ]
            },
        )


    def _write_pose_recovery_summary(
        self,
        *,
        pose_x: float,
    ) -> Path:
        poses_path = (
            self.root
            / "selected_poses.sdf"
        )

        molecule = Chem.MolFromSmiles("C")

        if molecule is None:
            raise RuntimeError(
                "Could not build test molecule."
            )

        conformer = Chem.Conformer(
            molecule.GetNumAtoms()
        )
        conformer.SetAtomPosition(
            0,
            Point3D(
                float(pose_x),
                0.0,
                0.0,
            ),
        )
        molecule.AddConformer(
            conformer,
            assignId=True,
        )

        writer = Chem.SDWriter(
            str(poses_path)
        )

        try:
            writer.write(molecule)
        finally:
            writer.close()

        return self._write_json(
            "pose_set_recovery_summary.json",
            {
                "poses_sdf": str(
                    poses_path
                ),
                "top_cnn_pose": {
                    "pose_index": 1,
                    "receptor_conformer_id": (
                        "submitted_receptor"
                    ),
                    "pocket_id": (
                        "fpocket_01_pocket_1"
                    ),
                    "seed": 123,
                    "source_pose_number": 1,
                },
            },
        )

    def _evaluate(
        self,
        *,
        receptor_x: float,
        include_outlier: bool = True,
        goals_met: bool = False,
        pose_x: float | None = None,
    ) -> dict:
        return evaluate_structure_pocket_quality(
            structure_path=(
                self._write_receptor(
                    x=receptor_x
                )
            ),
            ramachandran_report_path=(
                self._write_ramachandran_report(
                    include_outlier=(
                        include_outlier
                    ),
                    goals_met=goals_met,
                )
            ),
            pocket_definitions_path=(
                self._write_pocket_definitions()
            ),
            pocket_selection_summary_path=(
                self._write_selection_summary()
            ),
            pose_recovery_summary_path=(
                self._write_pose_recovery_summary(
                    pose_x=pose_x
                )
                if pose_x is not None
                else None
            ),
        )

    def test_outlier_inside_selected_box(
        self,
    ) -> None:
        report = self._evaluate(
            receptor_x=0.0
        )

        self.assertEqual(
            report["verdict"],
            "selected_pocket_geometry_concern",
        )
        self.assertEqual(
            report[
                "inside_selected_box_outliers"
            ],
            ["ALA:A:1"],
        )
        self.assertEqual(
            report[
                "selected_box_local_outliers"
            ],
            ["ALA:A:1"],
        )

        nearest = report["outliers"][0][
            "nearest_selected_pocket"
        ]

        self.assertIsNotNone(nearest)
        self.assertEqual(
            nearest["localization"],
            "inside_docking_box",
        )
        self.assertAlmostEqual(
            nearest[
                "minimum_distance_to_box_angstrom"
            ],
            0.0,
        )

    def test_outlier_near_selected_box(
        self,
    ) -> None:
        # The box extends from x=-10 to x=10.
        # A residue at x=12 is 2 Å outside.
        report = self._evaluate(
            receptor_x=12.0
        )

        self.assertEqual(
            report["verdict"],
            "manual_review_of_selected_pocket",
        )
        self.assertEqual(
            report[
                "near_selected_box_outliers"
            ],
            ["ALA:A:1"],
        )

        nearest = report["outliers"][0][
            "nearest_selected_pocket"
        ]

        self.assertEqual(
            nearest["localization"],
            "near_docking_box",
        )
        self.assertAlmostEqual(
            nearest[
                "minimum_distance_to_box_angstrom"
            ],
            2.0,
            places=6,
        )


    def test_near_box_remote_from_pose_is_box_edge_only(
        self,
    ) -> None:
        report = self._evaluate(
            receptor_x=12.0,
            pose_x=0.0,
        )

        self.assertEqual(
            report["verdict"],
            (
                "usable_with_global_"
                "geometry_caution"
            ),
        )
        self.assertEqual(
            report[
                "near_selected_box_outliers"
            ],
            ["ALA:A:1"],
        )
        self.assertEqual(
            report[
                "selected_box_local_outliers"
            ],
            ["ALA:A:1"],
        )
        self.assertEqual(
            report[
                "selected_pose_local_outliers"
            ],
            [],
        )
        self.assertEqual(
            report[
                "box_edge_only_outliers"
            ],
            ["ALA:A:1"],
        )

        outlier = report["outliers"][0]

        self.assertEqual(
            outlier[
                "selected_pose_localization"
            ],
            "box_edge_only",
        )
        self.assertAlmostEqual(
            outlier[
                "minimum_distance_to_"
                "selected_pose_angstrom"
            ],
            12.0,
            places=6,
        )

    def test_direct_pose_contact_is_geometry_concern(
        self,
    ) -> None:
        report = self._evaluate(
            receptor_x=12.0,
            pose_x=12.5,
        )

        self.assertEqual(
            report["verdict"],
            (
                "selected_pocket_geometry_"
                "concern"
            ),
        )
        self.assertEqual(
            report[
                "direct_selected_pose_"
                "contact_outliers"
            ],
            ["ALA:A:1"],
        )
        self.assertEqual(
            report[
                "selected_pose_local_outliers"
            ],
            ["ALA:A:1"],
        )
        self.assertEqual(
            report[
                "box_edge_only_outliers"
            ],
            [],
        )

        outlier = report["outliers"][0]

        self.assertEqual(
            outlier[
                "selected_pose_localization"
            ],
            "direct_selected_pose_contact",
        )
        self.assertAlmostEqual(
            outlier[
                "minimum_distance_to_"
                "selected_pose_angstrom"
            ],
            0.5,
            places=6,
        )

    def test_near_pose_requires_manual_review(
        self,
    ) -> None:
        report = self._evaluate(
            receptor_x=12.0,
            pose_x=17.0,
        )

        self.assertEqual(
            report["verdict"],
            "manual_review_of_selected_pocket",
        )
        self.assertEqual(
            report[
                "near_selected_pose_outliers"
            ],
            ["ALA:A:1"],
        )
        self.assertEqual(
            report[
                "selected_pose_local_outliers"
            ],
            ["ALA:A:1"],
        )

        outlier = report["outliers"][0]

        self.assertEqual(
            outlier[
                "selected_pose_localization"
            ],
            "near_selected_pose",
        )
        self.assertAlmostEqual(
            outlier[
                "minimum_distance_to_"
                "selected_pose_angstrom"
            ],
            5.0,
            places=6,
        )

    def test_distal_outlier_produces_global_caution(
        self,
    ) -> None:
        # A residue at x=25 is 15 Å beyond
        # the nearest box boundary.
        report = self._evaluate(
            receptor_x=25.0
        )

        self.assertEqual(
            report["verdict"],
            (
                "usable_with_global_"
                "geometry_caution"
            ),
        )
        self.assertEqual(
            report[
                "selected_box_local_outliers"
            ],
            [],
        )

        nearest = report["outliers"][0][
            "nearest_selected_pocket"
        ]

        self.assertEqual(
            nearest["localization"],
            "distal_from_docking_box",
        )
        self.assertAlmostEqual(
            nearest[
                "minimum_distance_to_box_angstrom"
            ],
            15.0,
            places=6,
        )

    def test_clean_global_geometry_is_strong(
        self,
    ) -> None:
        report = self._evaluate(
            receptor_x=0.0,
            include_outlier=False,
            goals_met=True,
        )

        self.assertEqual(
            report["verdict"],
            "strong",
        )
        self.assertTrue(
            report["global_goals_met"]
        )
        self.assertEqual(
            report["outlier_count"],
            0,
        )

    def test_runner_writes_report(
        self,
    ) -> None:
        receptor = self._write_receptor(
            x=25.0
        )
        ramachandran = (
            self._write_ramachandran_report()
        )
        pockets = (
            self._write_pocket_definitions()
        )
        selection = (
            self._write_selection_summary()
        )

        output_dir = self.root / "output"

        result = run_structure_pocket_quality(
            structure_path=receptor,
            ramachandran_report_path=(
                ramachandran
            ),
            pocket_definitions_path=pockets,
            pocket_selection_summary_path=(
                selection
            ),
            output_dir=output_dir,
        )

        output_path = (
            output_dir
            / "structure_pocket_quality.json"
        )

        self.assertEqual(
            result["output_path"],
            output_path,
        )
        self.assertTrue(
            output_path.is_file()
        )

        saved_report = json.loads(
            output_path.read_text(
                encoding="utf-8"
            )
        )

        self.assertEqual(
            saved_report["schema_version"],
            "structure_pocket_quality.v0.2",
        )
        self.assertEqual(
            saved_report["verdict"],
            (
                "usable_with_global_"
                "geometry_caution"
            ),
        )


if __name__ == "__main__":
    unittest.main()
