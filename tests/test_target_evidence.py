from __future__ import annotations

import unittest

from compoundrank.target_evidence import build_target_evidence


class TargetEvidenceTests(unittest.TestCase):
    def test_hiv_like_protease_inference(self) -> None:
        summary = {
            "job_id": "test",
            "status": "complete",
            "result_counts": {
                "cdd": 1,
                "interpro": 1,
                "vogdb": 0,
            },
            "rows": {
                "cdd": [
                    {
                        "name": "HIV-1 protease",
                        "description": "retropepsin aspartyl protease peptidase",
                    }
                ],
                "interpro": [
                    {
                        "signature_desc": "Aspartic peptidase family",
                    }
                ],
                "vogdb": [],
            },
        }

        evidence = build_target_evidence(summary)

        interpretation = evidence["target_interpretation"]
        future_query = evidence["future_ligand_database_query"]

        self.assertEqual(interpretation["target_class"], "viral protease")
        self.assertEqual(interpretation["enzyme_class"], "aspartyl protease")
        self.assertEqual(interpretation["docking_priority"], "high")
        self.assertEqual(interpretation["evidence_confidence"], "high")
        self.assertIn("retroviral", interpretation["target_name"].lower())
        self.assertIn(
            "HIV protease inhibitor",
            future_query["query_terms"],
        )

    def test_unknown_when_no_hits(self) -> None:
        summary = {
            "job_id": "empty",
            "status": "complete",
            "result_counts": {
                "cdd": 0,
                "interpro": 0,
                "vogdb": 0,
            },
            "rows": {
                "cdd": [],
                "interpro": [],
                "vogdb": [],
            },
        }

        evidence = build_target_evidence(summary)

        self.assertEqual(
            evidence["target_interpretation"]["target_class"],
            "unknown",
        )
        self.assertEqual(
            evidence["target_interpretation"]["docking_priority"],
            "low",
        )


if __name__ == "__main__":
    unittest.main()
