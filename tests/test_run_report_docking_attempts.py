import csv
import tempfile
import unittest
from pathlib import Path

from compoundrank.run_report import write_run_report


class RunReportDockingAttemptTests(unittest.TestCase):
    def test_run_report_includes_failed_docking_attempts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)

            with (output_dir / "docking_attempt_summary.csv").open(
                "w",
                newline="",
                encoding="utf-8",
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "compound",
                        "pocket",
                        "raw_poses",
                        "accepted_poses",
                        "rejected_poses",
                        "status",
                        "best_raw_cnn_score",
                        "best_accepted_cnn_score",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "compound": "darunavir",
                        "pocket": "user_box_01",
                        "raw_poses": "20",
                        "accepted_poses": "1",
                        "rejected_poses": "19",
                        "status": "accepted",
                        "best_raw_cnn_score": "0.815858543",
                        "best_accepted_cnn_score": "0.815858543",
                    }
                )
                writer.writerow(
                    {
                        "compound": "ritonavir",
                        "pocket": "user_box_01",
                        "raw_poses": "20",
                        "accepted_poses": "0",
                        "rejected_poses": "20",
                        "status": "failed_posebusters",
                        "best_raw_cnn_score": "0.624900000",
                        "best_accepted_cnn_score": "",
                    }
                )

            report_path = write_run_report(output_dir=output_dir)
            text = report_path.read_text(encoding="utf-8")

            self.assertIn("## Docking Attempt Summary", text)
            self.assertIn("### Attempted Ligands With No Accepted Final Pose", text)
            self.assertIn("ritonavir", text)
            self.assertIn("failed_posebusters", text)
            self.assertIn("0.624900000", text)


if __name__ == "__main__":
    unittest.main()
