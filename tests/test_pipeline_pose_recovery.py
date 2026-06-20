import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D

from compoundrank.models import PoseRecord
from compoundrank.pipeline import (
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

            second = molecule_at_offset(3.0)
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


if __name__ == "__main__":
    unittest.main()
