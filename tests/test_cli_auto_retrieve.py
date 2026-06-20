import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from compoundrank.cli import build_parser, main


class CliAutoRetrieveTests(unittest.TestCase):
    def test_parser_accepts_auto_retrieve_flags(self):
        args = build_parser().parse_args(
            [
                "--receptor", "/tmp/receptor.pdb",
                "--data-root", "/tmp/data",
                "--fasta", "/tmp/protein.faa",
                "--auto-retrieve-ligands",
                "--auto-retrieve-mode", "generic-strict",
                "--auto-retrieve-max-candidates", "4",
                "--auto-retrieve-no-fetch-structures",
            ]
        )

        self.assertTrue(args.auto_retrieve_ligands)
        self.assertEqual(args.auto_retrieve_mode, "generic-strict")
        self.assertEqual(args.auto_retrieve_max_candidates, 4)
        self.assertTrue(args.auto_retrieve_no_fetch_structures)

    def test_auto_retrieve_requires_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            receptor = tmp / "receptor.pdb"
            receptor.write_text("END\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                main(
                    [
                        "--receptor", str(receptor),
                        "--data-root", str(tmp / "data"),
                        "--auto-retrieve-ligands",
                    ]
                )

    def test_auto_retrieve_passes_empty_manual_requests_to_pipeline(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            receptor = tmp / "receptor.pdb"
            fasta = tmp / "protein.faa"
            receptor.write_text("END\n", encoding="utf-8")
            fasta.write_text(">protein\nAAAA\n", encoding="utf-8")

            with patch("compoundrank.cli.run_pipeline") as mocked:
                main(
                    [
                        "--receptor", str(receptor),
                        "--data-root", str(tmp / "data"),
                        "--fasta", str(fasta),
                        "--auto-retrieve-ligands",
                        "--auto-retrieve-mode", "generic-strict",
                        "--auto-retrieve-max-candidates", "3",
                    ]
                )

            kwargs = mocked.call_args.kwargs
            self.assertEqual(kwargs["ligand_requests"], [])
            self.assertTrue(kwargs["auto_retrieve_ligands"])
            self.assertEqual(
                kwargs["auto_retrieve_mode"],
                "generic-strict",
            )
            self.assertEqual(kwargs["auto_retrieve_max_candidates"], 3)


    def test_parser_accepts_pose_recovery_flags(
        self,
    ):
        args = build_parser().parse_args(
            [
                "--receptor",
                "/tmp/receptor.pdb",
                "--data-root",
                "/tmp/data",
                "--ligand-file",
                "/tmp/ligand.sdf",
                "--reference-ligand",
                "/tmp/reference.sdf",
                "--pose-recovery-rmsd-threshold",
                "1.75",
            ]
        )

        self.assertEqual(
            args.reference_ligand,
            "/tmp/reference.sdf",
        )
        self.assertAlmostEqual(
            args.pose_recovery_rmsd_threshold,
            1.75,
        )

    def test_pose_recovery_arguments_are_forwarded_to_pipeline(
        self,
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)

            receptor = tmp / "receptor.pdb"
            ligand = tmp / "ligand.sdf"
            reference = tmp / "reference.sdf"

            receptor.write_text(
                "END\n",
                encoding="utf-8",
            )
            ligand.write_text(
                "ligand\n",
                encoding="utf-8",
            )
            reference.write_text(
                "reference\n",
                encoding="utf-8",
            )

            with patch(
                "compoundrank.cli.run_pipeline"
            ) as mocked:
                main(
                    [
                        "--receptor",
                        str(receptor),
                        "--data-root",
                        str(tmp / "data"),
                        "--ligand-file",
                        str(ligand),
                        "--reference-ligand",
                        str(reference),
                        "--pose-recovery-rmsd-threshold",
                        "1.5",
                    ]
                )

            kwargs = mocked.call_args.kwargs

            self.assertEqual(
                kwargs["reference_ligand"],
                reference.resolve(),
            )
            self.assertAlmostEqual(
                kwargs[
                    "pose_recovery_rmsd_threshold"
                ],
                1.5,
            )


if __name__ == "__main__":
    unittest.main()
