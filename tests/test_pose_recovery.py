import json
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D

from compoundrank.pose_recovery import (
    compare_pose_recovery,
    evaluate_scored_pose_sdf,
    symmetry_aware_nofit_rmsd,
    write_outputs,
    write_scored_pose_outputs,
)


def pdb_line(record, serial, atom, resname, chain, resseq, x, y, z, element="C"):
    return (
        f"{record:<6}{serial:>5} "
        f"{atom:<4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"  1.00 20.00          {element:>2}"
    )


class PoseRecoveryTests(unittest.TestCase):
    def test_compare_pose_recovery_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)

            reference = tmp / "reference_ligand.pdb"
            reference.write_text(
                "\n".join(
                    [
                        pdb_line("HETATM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, 0.0),
                        pdb_line("HETATM", 2, "C2", "LIG", "A", 1, 2.0, 0.0, 0.0),
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            docked = tmp / "docked_pose.pdb"
            docked.write_text(
                "\n".join(
                    [
                        pdb_line("ATOM", 1, "CA", "PRO", "A", 1, -10.0, 0.0, 0.0),
                        pdb_line("HETATM", 2, "C1", "UNL", "A", 2, 1.0, 0.0, 0.0),
                        pdb_line("HETATM", 3, "C2", "UNL", "A", 2, 3.0, 0.0, 0.0),
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            metrics = compare_pose_recovery(reference_ligand=reference, docked_pose=docked, use_openbabel=False)

            self.assertEqual(metrics.reference_atom_count, 2)
            self.assertEqual(metrics.docked_atom_count, 2)
            self.assertAlmostEqual(metrics.center_distance, 1.0)
            self.assertAlmostEqual(metrics.ordered_coordinate_rmsd, 1.0)
            self.assertEqual(metrics.interpretation, "strong_pose_recovery")
            self.assertIsNone(metrics.openbabel_rmsd)
            self.assertIsNone(metrics.openbabel_minimized_rmsd)

            outputs = write_outputs(metrics, tmp / "out")
            self.assertTrue(outputs["pose_recovery_metrics_json"].exists())
            self.assertTrue(outputs["pose_recovery_metrics_csv"].exists())
            self.assertTrue(outputs["pose_recovery_report"].exists())

            data = json.loads(outputs["pose_recovery_metrics_json"].read_text(encoding="utf-8"))
            self.assertAlmostEqual(data["center_distance"], 1.0)


    def test_scored_sdf_pose_recovery_tolerates_bond_order_encoding(
        self,
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)

            def build_molecule(
                *,
                bond_type,
                translation,
                cnnscore=None,
            ):
                editable = Chem.RWMol()

                first = editable.AddAtom(
                    Chem.Atom("C")
                )
                second = editable.AddAtom(
                    Chem.Atom("C")
                )

                editable.AddBond(
                    first,
                    second,
                    bond_type,
                )

                molecule = editable.GetMol()

                conformer = Chem.Conformer(2)
                conformer.Set3D(True)
                conformer.SetAtomPosition(
                    0,
                    Point3D(
                        0.0 + translation,
                        0.0,
                        0.0,
                    ),
                )
                conformer.SetAtomPosition(
                    1,
                    Point3D(
                        1.5 + translation,
                        0.0,
                        0.0,
                    ),
                )

                molecule.AddConformer(
                    conformer,
                    assignId=True,
                )

                if cnnscore is not None:
                    molecule.SetProp(
                        "CNNscore",
                        str(cnnscore),
                    )
                    molecule.SetProp(
                        "CNNaffinity",
                        "7.0",
                    )
                    molecule.SetProp(
                        "minimizedAffinity",
                        "-8.0",
                    )

                return molecule

            reference = build_molecule(
                bond_type=(
                    Chem.BondType.SINGLE
                ),
                translation=0.0,
            )

            top_pose = build_molecule(
                bond_type=(
                    Chem.BondType.DOUBLE
                ),
                translation=1.0,
                cnnscore=0.9,
            )

            lower_pose = build_molecule(
                bond_type=(
                    Chem.BondType.DOUBLE
                ),
                translation=3.0,
                cnnscore=0.5,
            )

            direct_rmsd, mapping_count = (
                symmetry_aware_nofit_rmsd(
                    top_pose,
                    reference,
                )
            )

            self.assertAlmostEqual(
                direct_rmsd,
                1.0,
            )
            self.assertGreaterEqual(
                mapping_count,
                1,
            )

            reference_path = (
                tmp / "reference.sdf"
            )
            poses_path = (
                tmp / "poses.sdf"
            )

            reference_writer = (
                Chem.SDWriter(
                    str(reference_path)
                )
            )
            reference_writer.write(
                reference
            )
            reference_writer.close()

            poses_writer = Chem.SDWriter(
                str(poses_path)
            )
            poses_writer.write(
                top_pose
            )
            poses_writer.write(
                lower_pose
            )
            poses_writer.close()

            summary = (
                evaluate_scored_pose_sdf(
                    reference_ligand=(
                        reference_path
                    ),
                    poses_sdf=poses_path,
                    rmsd_threshold=2.0,
                )
            )

            self.assertEqual(
                summary["mapped_pose_count"],
                2,
            )
            self.assertEqual(
                summary["top_cnn_pose"][
                    "pose_index"
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
            self.assertEqual(
                summary["overall"],
                (
                    "cognate_pose_recovery_"
                    "and_ranking_pass"
                ),
            )

            outputs = (
                write_scored_pose_outputs(
                    summary,
                    tmp / "batch_out",
                )
            )

            for output_path in (
                outputs.values()
            ):
                self.assertTrue(
                    output_path.exists()
                )


if __name__ == "__main__":
    unittest.main()
