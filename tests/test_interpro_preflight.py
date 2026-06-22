import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cpu_server.tools import interpro_tool


class InterProDatabasePreflightTests(unittest.TestCase):
    def test_parse_applications_normalizes_and_deduplicates(self):
        self.assertEqual(
            interpro_tool._parse_applications(
                " PFAM,ncbifam,pfam, SUPERFAMILY "
            ),
            ["pfam", "ncbifam", "superfamily"],
        )

    def test_missing_database_fails_and_writes_report(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            datadir = root / "data"
            output_dir = root / "output"
            datadir.mkdir()

            with self.assertRaises(
                interpro_tool.InterProDatabasePreflightError
            ) as context:
                interpro_tool.preflight_interpro_databases(
                    datadir=datadir,
                    output_dir=output_dir,
                    applications="pfam,ncbifam",
                    mode="local",
                )

            message = str(context.exception)
            self.assertIn("database_missing", message)
            self.assertIn("pfam", message)
            self.assertIn("ncbifam", message)

            report_path = (
                output_dir
                / "interpro_database_preflight.json"
            )
            self.assertTrue(report_path.exists())

            report = json.loads(
                report_path.read_text(encoding="utf-8")
            )

            self.assertEqual(report["status"], "failed")
            self.assertEqual(
                {
                    check["application"]
                    for check in report["checks"]
                },
                {"pfam", "ncbifam"},
            )
            self.assertTrue(
                all(
                    check["error_type"]
                    == "database_missing"
                    for check in report["checks"]
                )
            )

    @patch(
        "cpu_server.tools.interpro_tool.shutil.which",
        return_value="/usr/bin/hmmstat",
    )
    @patch(
        "cpu_server.tools.interpro_tool.subprocess.run"
    )
    def test_invalid_hmm_is_rejected(
        self,
        mock_run,
        mock_which,
    ):
        mock_run.return_value = subprocess.CompletedProcess(
            args=["hmmstat"],
            returncode=1,
            stdout="",
            stderr=(
                "Error: bad file format in HMM file "
                "19.0/ncbifam.hmm"
            ),
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            datadir = root / "data"
            output_dir = root / "output"
            hmm_path = (
                datadir
                / "19.0"
                / "ncbifam.hmm"
            )
            hmm_path.parent.mkdir(parents=True)
            hmm_path.write_text(
                "not a valid HMM\n",
                encoding="utf-8",
            )

            with self.assertRaises(
                interpro_tool.InterProDatabasePreflightError
            ) as context:
                interpro_tool.preflight_interpro_databases(
                    datadir=datadir,
                    output_dir=output_dir,
                    applications="ncbifam",
                    mode="local",
                )

            message = str(context.exception)
            self.assertIn("database_invalid", message)
            self.assertIn("bad file format", message)

            report = json.loads(
                (
                    output_dir
                    / "interpro_database_preflight.json"
                ).read_text(encoding="utf-8")
            )

            check = report["checks"][0]
            self.assertEqual(check["status"], "failed")
            self.assertEqual(
                check["error_type"],
                "database_invalid",
            )
            self.assertEqual(
                check["hmmstat"]["return_code"],
                1,
            )
            self.assertIn(
                "bad file format",
                check["hmmstat"]["stderr"],
            )

            mock_which.assert_called_once_with("hmmstat")
            mock_run.assert_called_once()

    @patch(
        "cpu_server.tools.interpro_tool.shutil.which",
        return_value="/usr/bin/hmmstat",
    )
    @patch(
        "cpu_server.tools.interpro_tool.subprocess.run"
    )
    def test_valid_hmm_passes_and_records_selected_file(
        self,
        mock_run,
        mock_which,
    ):
        mock_run.return_value = subprocess.CompletedProcess(
            args=["hmmstat"],
            returncode=0,
            stdout="# idx name accession\n1 model -\n",
            stderr="",
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            datadir = root / "data"
            output_dir = root / "output"

            older = (
                datadir
                / "18.0"
                / "ncbifam.hmm"
            )
            newer = (
                datadir
                / "19.0"
                / "ncbifam.hmm"
            )

            older.parent.mkdir(parents=True)
            newer.parent.mkdir(parents=True)

            older.write_text(
                "HMMER3/f\nolder\n",
                encoding="utf-8",
            )
            newer.write_text(
                "HMMER3/f\nnewer\n",
                encoding="utf-8",
            )

            report = (
                interpro_tool
                .preflight_interpro_databases(
                    datadir=datadir,
                    output_dir=output_dir,
                    applications="ncbifam",
                    mode="local",
                )
            )

            self.assertEqual(report["status"], "passed")
            self.assertEqual(
                report["checks"][0]["selected_path"],
                str(newer.resolve()),
            )
            self.assertEqual(
                report["checks"][0]["status"],
                "passed",
            )

            mock_which.assert_called_once_with("hmmstat")
            mock_run.assert_called_once()

    @patch(
        "cpu_server.tools.interpro_tool._run"
    )
    def test_failed_preflight_prevents_nextflow_launch(
        self,
        mock_nextflow_run,
    ):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            fasta_path = root / "protein.fasta"
            output_dir = root / "output"
            datadir = root / "data"

            fasta_path.write_text(
                ">test\nMKTIIALSYIFCLVFA\n",
                encoding="utf-8",
            )
            datadir.mkdir()

            with self.assertRaises(
                interpro_tool.InterProDatabasePreflightError
            ):
                interpro_tool.run_interpro_local(
                    fasta_path=fasta_path,
                    output_dir=output_dir,
                    datadir=datadir,
                    applications="ncbifam",
                )

            mock_nextflow_run.assert_not_called()

            self.assertTrue(
                (
                    output_dir
                    / "interpro_database_preflight.json"
                ).exists()
            )


if __name__ == "__main__":
    unittest.main()
