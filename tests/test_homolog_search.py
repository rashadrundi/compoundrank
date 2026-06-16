from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from compoundrank.homolog_search import parse_cpu_response, run_homolog_search


class HomologSearchTests(unittest.TestCase):
    def test_parse_cpu_response_counts_rows(self) -> None:
        response = {
            "job_id": "abc",
            "status": "complete",
            "results": {
                "cdd": [{"id": 1}, {"id": 2}],
                "interpro": [{"id": 3}],
                "vogdb": [],
            },
            "files": {
                "report": "example.json",
            },
        }

        parsed = parse_cpu_response(response)

        self.assertEqual(parsed["job_id"], "abc")
        self.assertEqual(parsed["status"], "complete")
        self.assertEqual(parsed["result_counts"]["cdd"], 2)
        self.assertEqual(parsed["result_counts"]["interpro"], 1)
        self.assertEqual(parsed["result_counts"]["vogdb"], 0)
        self.assertEqual(parsed["files"]["report"], "example.json")


    def test_run_homolog_search_writes_target_evidence_outputs(self) -> None:
        fake_response = {
            "job_id": "test-job",
            "status": "complete",
            "result_counts": {
                "cdd": 1,
                "interpro": 0,
                "vogdb": 0,
            },
            "rows": {
                "cdd": [
                    {
                        "source": "CDD",
                        "name": "RVP",
                        "accession": "pfam00077",
                        "start": "5",
                        "end": "98",
                        "evalue": "5.0e-36",
                        "score": "117.4",
                        "notes": "specific retroviral protease domain",
                    }
                ],
                "interpro": [],
                "vogdb": [],
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "protein.faa"
            output_dir = root / "out"

            fasta.write_text(
                ">protein\n"
                "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF\n",
                encoding="utf-8",
            )

            with patch(
                "compoundrank.homolog_search.post_fasta",
                return_value=fake_response,
            ):
                result = run_homolog_search(
                    api_url="http://example.invalid/analyze/fasta",
                    fasta_path=fasta,
                    output_dir=output_dir,
                    timeout_seconds=1,
                )

            self.assertEqual(result["status"], "ok")
            self.assertTrue((output_dir / "homolog_search_raw.json").exists())
            self.assertTrue((output_dir / "homolog_search_summary.json").exists())
            self.assertTrue((output_dir / "target_evidence.json").exists())
            self.assertTrue((output_dir / "target_evidence_report.md").exists())

            evidence = json.loads(
                (output_dir / "target_evidence.json").read_text(
                    encoding="utf-8",
                )
            )

            self.assertEqual(
                evidence["target_interpretation"]["target_class"],
                "viral protease",
            )
            self.assertEqual(
                evidence["target_interpretation"]["evidence_confidence"],
                "high",
            )
            self.assertIn(
                "pfam00077",
                evidence["evidence"]["matched_keywords"],
            )


if __name__ == "__main__":
    unittest.main()
