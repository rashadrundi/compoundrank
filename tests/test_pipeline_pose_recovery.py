import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D

from compoundrank.models import PoseRecord
from compoundrank.pipeline import (
    _accepted_records_for_selected_pocket,
    _raw_records_for_selected_pocket_across_conformers,
    _run_selected_pocket_pose_recovery,
)


def molecule_at_offset(
    offset: float,
) -> Chem.Mol:
    molecule = Chem.MolFromSmiles("CC")

    if molecule is None:
        raise RuntimeError(
            "Could not create the test molecule."
        )

    conformer = Chem.Conformer(
        molecule.GetNumAtoms()
    )
    conformer.Set3D(True)

    conformer.SetAtomPosition(
        0,
        Point3D(offset, 0.0, 0.0),
    )
    conformer.SetAtomPosition(
        1,
        Point3D(offset + 1.5, 0.0, 0.0),
    )

    molecule.RemoveAllConformers()
    molecule.AddConformer(
        conformer,
        assignId=True,
    )

    return molecule


class PipelinePoseRecoveryTests(
    unittest.TestCase
):
    def test_selected_pocket_records_are_persisted_and_evaluated(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            reference_path = (
                root / "reference.sdf"
            )

            reference_writer = Chem.SDWriter(
                str(reference_path)
            )
            reference_writer.write(
                molecule_at_offset(0.0)
            )
            reference_writer.close()

            first = molecule_at_offset(1.0)
            first.SetDoubleProp(
                "CNNscore",
                0.9,
            )
            first.SetDoubleProp(
                "CNNaffinity",
                7.0,
            )
            first.SetDoubleProp(
                "minimizedAffinity",
                -8.0,
            )

            second = molecule_at_offset(0.5)
            second.SetDoubleProp(
                "CNNscore",
                0.5,
            )
            second.SetDoubleProp(
                "CNNaffinity",
                6.0,
            )
            second.SetDoubleProp(
                "minimizedAffinity",
                -7.0,
            )

            records = [
                PoseRecord(
                    ligand_name="example",
                    seed=100,
                    pose_number=1,
                    molecule=first,
                    cnn_score=0.9,
                    cnn_affinity=7.0,
                    minimized_affinity=-8.0,
                    source_sdf=Path(
                        "seed_100/poses.sdf"
                    ),
                    pocket_id="fpocket_02",
                    pocket_rank=2,
                    fpocket_score=0.2,
                    receptor_conformer_id=(
                        "snapshot_0001"
                    ),
                ),
                PoseRecord(
                    ligand_name="example",
                    seed=200,
                    pose_number=1,
                    molecule=second,
                    cnn_score=0.5,
                    cnn_affinity=6.0,
                    minimized_affinity=-7.0,
                    source_sdf=Path(
                        "seed_200/poses.sdf"
                    ),
                    pocket_id="fpocket_02",
                    pocket_rank=2,
                    fpocket_score=0.2,
                    receptor_conformer_id=(
                        "submitted_receptor"
                    ),
                ),
            ]

            summary, outputs = (
                _run_selected_pocket_pose_recovery(
                    reference_ligand=(
                        reference_path
                    ),
                    records=records,
                    output_dir=root,
                    ligand_name="example",
                    pocket_id="fpocket_02",
                    rmsd_threshold=2.0,
                    selected_receptor_conformer_id=(
                        "snapshot_0001"
                    ),
                )
            )

            self.assertEqual(
                summary["mapped_pose_count"],
                2,
            )
            self.assertEqual(
                summary["evaluated_compound"],
                "example",
            )
            self.assertEqual(
                summary["evaluated_pocket_id"],
                "fpocket_02",
            )
            self.assertFalse(
                summary[
                    "reference_ligand_used_for_pocket_selection"
                ]
            )
            self.assertFalse(
                summary[
                    "reference_ligand_used_for_docking"
                ]
            )
            self.assertEqual(
                summary[
                    "normally_selected_receptor_conformer_id"
                ],
                "snapshot_0001",
            )
            self.assertEqual(
                summary[
                    "evaluated_receptor_conformer_ids"
                ],
                [
                    "snapshot_0001",
                    "submitted_receptor",
                ],
            )
            self.assertEqual(
                summary[
                    "evaluated_receptor_conformer_count"
                ],
                2,
            )
            self.assertEqual(
                summary["top_cnn_pose"][
                    "receptor_conformer_id"
                ],
                "snapshot_0001",
            )
            self.assertEqual(
                summary["best_sampled_pose"][
                    "receptor_conformer_id"
                ],
                "submitted_receptor",
            )
            self.assertAlmostEqual(
                summary["best_sampled_pose"][
                    "heavy_atom_rmsd"
                ],
                0.5,
            )
            self.assertEqual(
                summary["top_cnn_pose"]["seed"],
                100,
            )
            self.assertEqual(
                summary["top_cnn_pose"][
                    "source_pose_number"
                ],
                1,
            )
            self.assertAlmostEqual(
                summary["top_cnn_pose"][
                    "heavy_atom_rmsd"
                ],
                1.0,
            )
            self.assertTrue(
                summary["sampling_pass"]
            )
            self.assertTrue(
                summary["ranking_pass"]
            )

            persisted_pose_set = (
                root
                / "pose_recovery_selected_pocket_poses.sdf"
            )

            self.assertTrue(
                persisted_pose_set.exists()
            )
            self.assertGreater(
                persisted_pose_set.stat().st_size,
                0,
            )

            for output_path in outputs.values():
                self.assertTrue(
                    output_path.exists()
                )

    def test_only_selected_pocket_records_feed_hypotheses(
        self,
    ) -> None:
        selected_record = PoseRecord(
            ligand_name="example",
            seed=100,
            pose_number=1,
            molecule=None,
            cnn_score=0.30,
            cnn_affinity=4.0,
            minimized_affinity=-5.0,
            source_sdf=Path(
                "selected/poses.sdf"
            ),
            pocket_id="selected_pocket",
        )

        false_high_score_record = PoseRecord(
            ligand_name="example",
            seed=100,
            pose_number=1,
            molecule=None,
            cnn_score=0.99,
            cnn_affinity=9.0,
            minimized_affinity=-8.0,
            source_sdf=Path(
                "false/poses.sdf"
            ),
            pocket_id="false_pocket",
        )

        records = (
            _accepted_records_for_selected_pocket(
                valid_records_by_pocket={
                    "selected_pocket": [
                        selected_record
                    ],
                    "false_pocket": [
                        false_high_score_record
                    ],
                },
                selected_pocket_id=(
                    "selected_pocket"
                ),
            )
        )

        self.assertEqual(
            records,
            [selected_record],
        )

        self.assertNotIn(
            false_high_score_record,
            records,
        )


    def test_raw_selected_pocket_records_include_all_conformers(
        self,
    ) -> None:
        submitted = PoseRecord(
            ligand_name="example",
            seed=100,
            pose_number=1,
            molecule=None,
            cnn_score=0.50,
            cnn_affinity=5.0,
            minimized_affinity=-6.0,
            source_sdf=Path("submitted/poses.sdf"),
            pocket_id="selected_pocket",
            receptor_conformer_id="submitted_receptor",
        )
        snapshot = PoseRecord(
            ligand_name="example",
            seed=200,
            pose_number=1,
            molecule=None,
            cnn_score=0.60,
            cnn_affinity=6.0,
            minimized_affinity=-7.0,
            source_sdf=Path("snapshot/poses.sdf"),
            pocket_id="selected_pocket",
            receptor_conformer_id="snapshot_0001",
        )
        decoy = PoseRecord(
            ligand_name="example",
            seed=200,
            pose_number=2,
            molecule=None,
            cnn_score=0.99,
            cnn_affinity=9.0,
            minimized_affinity=-9.0,
            source_sdf=Path("decoy/poses.sdf"),
            pocket_id="decoy_pocket",
            receptor_conformer_id="snapshot_0001",
        )

        records = (
            _raw_records_for_selected_pocket_across_conformers(
                raw_records_by_pocket={
                    (
                        "submitted_receptor",
                        "selected_pocket",
                    ): [submitted],
                    (
                        "snapshot_0001",
                        "selected_pocket",
                    ): [snapshot],
                    (
                        "snapshot_0001",
                        "decoy_pocket",
                    ): [decoy],
                },
                selected_pocket_id="selected_pocket",
            )
        )

        self.assertEqual(
            records,
            [submitted, snapshot],
        )
        self.assertNotIn(decoy, records)



if __name__ == "__main__":
    unittest.main()
