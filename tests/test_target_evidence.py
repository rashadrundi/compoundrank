from __future__ import annotations

import unittest

from compoundrank.target_evidence import (
    build_target_evidence,
    render_target_evidence_report,
)


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
                        "source": "CDD",
                        "name": "RVP",
                        "accession": "pfam00077",
                        "start": "5",
                        "end": "98",
                        "evalue": "5.77331e-36",
                        "score": "117.468",
                        "notes": (
                            "specific retroviral protease "
                            "domain; retropepsin aspartyl "
                            "protease"
                        ),
                    }
                ],
                "interpro": [
                    {
                        "signature_desc": (
                            "Aspartic peptidase family"
                        ),
                    }
                ],
                "vogdb": [],
            },
        }

        evidence = build_target_evidence(summary)

        interpretation = evidence[
            "target_interpretation"
        ]
        future_query = evidence[
            "future_ligand_database_query"
        ]

        self.assertEqual(
            interpretation["target_class"],
            "viral protease",
        )
        self.assertEqual(
            interpretation["enzyme_class"],
            "aspartyl protease",
        )
        self.assertEqual(
            interpretation["docking_priority"],
            "high",
        )
        self.assertEqual(
            interpretation["evidence_confidence"],
            "high",
        )
        self.assertIn(
            "retroviral",
            interpretation["target_name"].lower(),
        )
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
            evidence[
                "target_interpretation"
            ]["target_class"],
            "unknown",
        )
        self.assertEqual(
            evidence[
                "target_interpretation"
            ]["docking_priority"],
            "low",
        )

    def test_generic_protease_does_not_imply_aspartyl_or_hiv(
        self,
    ) -> None:
        summary = {
            "job_id": "generic-protease",
            "status": "complete",
            "result_counts": {
                "cdd": 0,
                "interpro": 1,
                "vogdb": 1,
            },
            "rows": {
                "cdd": [],
                "interpro": [
                    {
                        "signature_desc": (
                            "Putative viral protease family"
                        ),
                    }
                ],
                "vogdb": [
                    {
                        "description": (
                            "Unclassified viral peptidase"
                        ),
                    }
                ],
            },
        }

        evidence = build_target_evidence(
            summary
        )

        interpretation = evidence[
            "target_interpretation"
        ]
        query_terms = evidence[
            "future_ligand_database_query"
        ]["query_terms"]

        self.assertEqual(
            interpretation["target_class"],
            "viral protease",
        )
        self.assertIsNone(
            interpretation["enzyme_class"]
        )
        self.assertFalse(
            any(
                "hiv" in term.lower()
                for term in query_terms
            )
        )
        self.assertFalse(
            any(
                "aspartyl" in term.lower()
                for term in query_terms
            )
        )

    def test_coronavirus_main_protease_inference(
        self,
    ) -> None:
        summary = {
            "job_id": "mpro",
            "status": "complete",
            "result_counts": {
                "cdd": 1,
                "interpro": 1,
                "vogdb": 1,
            },
            "rows": {
                "cdd": [
                    {
                        "name": "3C-like protease",
                        "accession": "benchmark",
                        "start": "1",
                        "end": "306",
                        "evalue": "1e-80",
                        "score": "250",
                    }
                ],
                "interpro": [
                    {
                        "signature_desc": (
                            "Peptidase C30, "
                            "coronavirus main protease"
                        ),
                    }
                ],
                "vogdb": [
                    {
                        "description": (
                            "Coronavirus main proteinase"
                        ),
                    }
                ],
            },
        }

        evidence = build_target_evidence(
            summary
        )

        interpretation = evidence[
            "target_interpretation"
        ]
        query_terms = evidence[
            "future_ligand_database_query"
        ]["query_terms"]
        special = evidence["evidence"][
            "special_domain_evidence"
        ]

        self.assertEqual(
            interpretation["target_class"],
            "viral protease",
        )
        self.assertEqual(
            interpretation["enzyme_class"],
            "cysteine protease",
        )
        self.assertIn(
            "coronavirus",
            interpretation["target_name"].lower(),
        )
        self.assertIn(
            "coronaviridae",
            interpretation["viral_family"].lower(),
        )
        self.assertEqual(
            interpretation["evidence_confidence"],
            "high",
        )
        self.assertIsNotNone(special)
        self.assertIn(
            "3CLpro inhibitor",
            query_terms,
        )
        self.assertFalse(
            any(
                "hiv" in term.lower()
                for term in query_terms
            )
        )

    def test_failed_interpro_is_not_reported_as_no_hits(
        self,
    ) -> None:
        error = (
            "RuntimeError: [InterPro] Command failed\n"
            "Exit code: 1\n"
            "Error: bad file format in HMM file "
            "19.0/ncbifam.hmm"
        )

        summary = {
            "job_id": "borna-regression",
            "status": "partial",
            "result_counts": {
                "cdd": 1,
                "interpro": 0,
                "vogdb": 6,
            },
            "tool_statuses": {
                "cdd": "complete",
                "interpro": "failed",
                "vogdb": "complete",
            },
            "tool_errors": {
                "interpro": error,
            },
            "rows": {
                "cdd": [
                    {
                        "source": "CDD",
                        "name": "Mononeg_RNA_pol",
                        "accession": "cl15638",
                        "start": "172",
                        "end": "776",
                        "evalue": "1.41641e-60",
                        "score": "222.217",
                        "notes": (
                            "mononegavirus RNA polymerase"
                        ),
                    }
                ],
                "interpro": [],
                "vogdb": [
                    {
                        "description": (
                            "viral RNA polymerase protein"
                        ),
                    }
                    for _ in range(6)
                ],
            },
        }

        evidence = build_target_evidence(
            summary,
            source_fasta="BoDV1_test.fasta",
        )

        source = evidence["source"]

        self.assertEqual(
            source["status"],
            "partial",
        )
        self.assertEqual(
            source["tool_statuses"]["interpro"],
            "failed",
        )
        self.assertEqual(
            source["result_counts"]["interpro"],
            0,
        )
        self.assertIn(
            "bad file format",
            source["tool_errors"]["interpro"],
        )

        reasoning = "\n".join(
            evidence["evidence"][
                "confidence_reasoning"
            ]
        )

        self.assertIn(
            "InterPro failed",
            reasoning,
        )
        self.assertIn(
            "0 usable result(s) were retained",
            reasoning,
        )
        self.assertNotIn(
            "InterPro completed successfully "
            "and returned no hits",
            reasoning,
        )

        limitations = "\n".join(
            evidence["limitations"]
        )

        self.assertIn(
            "InterPro failed",
            limitations,
        )
        self.assertIn(
            "must not be interpreted as a "
            "successful no-hit result",
            limitations,
        )

        report = render_target_evidence_report(
            evidence
        )

        self.assertIn(
            "| InterPro | failed | 0 |",
            report,
        )
        self.assertIn(
            "bad file format in HMM file "
            "19.0/ncbifam.hmm",
            report,
        )
        self.assertNotIn(
            "| InterPro | complete_no_hits | 0 |",
            report,
        )


if __name__ == "__main__":
    unittest.main()
