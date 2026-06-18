import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.pose_recovery import compare_pose_recovery, write_outputs


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


if __name__ == "__main__":
    unittest.main()
