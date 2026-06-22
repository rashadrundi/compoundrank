from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from compoundrank.homolog_search import (
    parse_cpu_response,
    run_homolog_search,
)


class HomologSearchTests(unittest.TestCase):
    def test_parse_cpu_response_counts_rows(self) -> None:
        response = {
            "job_id": "abc",
            "status": "complete",
            "results": {
                "cdd": [
                    {"id": 1},
                    {"id": 2},
                ],
                "interpro": [
                    {"id": 3},
                ],
                "vogdb": [],
            },
            "files": {
                "report": "example.json",
            },
        }

        parsed = parse_cpu_response(response)

        self.assertEqual(
            parsed["job_id"],
            "abc",
        )
        self.assertEqual(
            parsed["status"],
            "complete",
        )
        self.assertEqual(
            parsed["source_status"],
            "complete",
        )

        self.assertEqual(
            parsed["result_counts"]["cdd"],
            2,
        )
        self.assertEqual(
            parsed["result_counts"]["interpro"],
            1,
        )
        self.assertEqual(
            parsed["result_counts"]["vogdb"],
            0,
        )

        self.assertEqual(
            parsed["tool_statuses"]["cdd"],
            "complete",
        )
        self.assertEqual(
            parsed["tool_statuses"]["interpro"],
            "complete",
        )
        self.assertEqual(
            parsed["tool_statuses"]["vogdb"],
            "complete_no_hits",
        )

        self.assertEqual(
            parsed["tool_errors"],
            {},
        )
        self.assertEqual(
            parsed["files"]["report"],
            "example.json",
        )

    def test_parse_cpu_response_preserves_failed_tool_status(
        self,
    ) -> None:
        interpro_error = (
            "RuntimeError: [InterPro] Command failed\n"
            "Exit code: 1\n"
            "Error: bad file format in HMM file "
            "19.0/ncbifam.hmm"
        )

        response = {
            "job_id": "borna-regression",
            "status": "partial_or_failed",
            "rows": {
                "cdd": [
                    {
                        "name": "Mononeg_RNA_pol",
                    }
                ],
                "interpro": [],
                "vogdb": [
                    {"id": index}
                    for index in range(6)
                ],
            },
            "cdd": {
                "status": "complete",
                "rows": [
                    {
                        "name": "Mononeg_RNA_pol",
                    }
                ],
            },
            "interpro": {
                "status": "failed",
                "rows": [],
                "error": interpro_error,
                "command": [
                    "run_interpro_tool",
                    "--applications",
                    "pfam,ncbifam",
                ],
            },
            "vogdb": {
                "status": "complete",
                "rows": [
                    {"id": index}
                    for index in range(6)
                ],
            },
        }

        parsed = parse_cpu_response(response)

        self.assertEqual(
            parsed["status"],
            "partial",
        )
        self.assertEqual(
            parsed["source_status"],
            "partial_or_failed",
        )

        self.assertEqual(
            parsed["result_counts"],
            {
                "cdd": 1,
                "interpro": 0,
                "vogdb": 6,
            },
        )

        self.assertEqual(
            parsed["tool_statuses"],
            {
                "cdd": "complete",
                "interpro": "failed",
                "vogdb": "complete",
            },
        )

        self.assertIn(
            "interpro",
            parsed["tool_errors"],
        )
        self.assertIn(
            "bad file format",
            parsed["tool_errors"]["interpro"],
        )

        self.assertEqual(
            parsed["tools"]["interpro"]["row_count"],
            0,
        )
        self.assertEqual(
            parsed["tools"]["interpro"]["status"],
            "failed",
        )
        self.assertEqual(
            parsed["tools"]["interpro"]["command"],
            [
                "run_interpro_tool",
                "--applications",
                "pfam,ncbifam",
            ],
        )

        # Critical regression assertion:
        # zero rows from a failed tool must never become
        # a successful no-hit result.
        self.assertNotEqual(
            parsed["tool_statuses"]["interpro"],
            "complete_no_hits",
        )

    def test_parse_cpu_response_distinguishes_successful_no_hits(
        self,
    ) -> None:
        response = {
            "job_id": "true-no-hits",
            "status": "complete",
            "rows": {
                "cdd": [],
                "interpro": [],
                "vogdb": [],
            },
            "cdd": {
                "status": "complete_no_hits",
                "rows": [],
            },
            "interpro": {
                "status": "complete_no_hits",
                "rows": [],
            },
            "vogdb": {
                "status": "complete_no_hits",
                "rows": [],
            },
        }

        parsed = parse_cpu_response(response)

        self.assertEqual(
            parsed["status"],
            "complete",
        )
        self.assertEqual(
            parsed["result_counts"],
            {
                "cdd": 0,
                "interpro": 0,
                "vogdb": 0,
            },
        )
        self.assertEqual(
            parsed["tool_statuses"],
            {
                "cdd": "complete_no_hits",
                "interpro": "complete_no_hits",
                "vogdb": "complete_no_hits",
            },
        )

    def test_unknown_status_is_not_converted_to_no_hits(
        self,
    ) -> None:
        response = {
            "job_id": "legacy-partial",
            "status": "partial_or_failed",
            "rows": {
                "cdd": [],
                "interpro": [],
                "vogdb": [],
            },
        }

        parsed = parse_cpu_response(response)

        self.assertEqual(
            parsed["status"],
            "unknown",
        )

        self.assertEqual(
            parsed["tool_statuses"],
            {
                "cdd": "unknown",
                "interpro": "unknown",
                "vogdb": "unknown",
            },
        )

        self.assertNotIn(
            "complete_no_hits",
            parsed["tool_statuses"].values(),
        )

    def test_run_homolog_search_writes_target_evidence_outputs(
        self,
    ) -> None:
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
                        "notes": (
                            "specific retroviral "
                            "protease domain"
                        ),
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
                "PQITLWQRPLVTIKIGGQLKEALLDTGADDTV"
                "LEEMNLPGRWKPKMIGGIGGFIKVRQYDQILI"
                "EICGHKAIGTVLVGPTPVNIIGRNLLTQIGCT"
                "LNF\n",
                encoding="utf-8",
            )

            with patch(
                "compoundrank.homolog_search.post_fasta",
                return_value=fake_response,
            ):
                result = run_homolog_search(
                    api_url=(
                        "http://example.invalid/"
                        "analyze/fasta"
                    ),
                    fasta_path=fasta,
                    output_dir=output_dir,
                    timeout_seconds=1,
                )

            self.assertEqual(
                result["status"],
                "ok",
            )
            self.assertEqual(
                result["cpu_status"],
                "complete",
            )
            self.assertEqual(
                result["tool_statuses"],
                {
                    "cdd": "complete",
                    "interpro": "complete_no_hits",
                    "vogdb": "complete_no_hits",
                },
            )

            self.assertTrue(
                (
                    output_dir
                    / "homolog_search_raw.json"
                ).exists()
            )
            self.assertTrue(
                (
                    output_dir
                    / "homolog_search_summary.json"
                ).exists()
            )
            self.assertTrue(
                (
                    output_dir
                    / "target_evidence.json"
                ).exists()
            )
            self.assertTrue(
                (
                    output_dir
                    / "target_evidence_report.md"
                ).exists()
            )

            summary = json.loads(
                (
                    output_dir
                    / "homolog_search_summary.json"
                ).read_text(
                    encoding="utf-8",
                )
            )

            self.assertEqual(
                summary["tool_statuses"]["cdd"],
                "complete",
            )
            self.assertEqual(
                summary["tool_statuses"]["interpro"],
                "complete_no_hits",
            )

            evidence = json.loads(
                (
                    output_dir
                    / "target_evidence.json"
                ).read_text(
                    encoding="utf-8",
                )
            )

            self.assertEqual(
                evidence[
                    "target_interpretation"
                ]["target_class"],
                "viral protease",
            )
            self.assertEqual(
                evidence[
                    "target_interpretation"
                ]["evidence_confidence"],
                "high",
            )
            self.assertIn(
                "pfam00077",
                evidence["evidence"][
                    "matched_keywords"
                ],
            )


if __name__ == "__main__":
    unittest.main()
