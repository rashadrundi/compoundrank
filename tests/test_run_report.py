from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.run_report import (
    _render_target_section,
    collect_pdb_hypotheses,
    write_run_report,
)


class RunReportTests(unittest.TestCase):
    def test_collect_pdb_hypotheses_from_filename_and_remarks(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            pdb = (
                root
                / (
                    "01__darunavir__"
                    "fpocket_02_pocket_2__"
                    "hypothesis_01.pdb"
                )
            )

            pdb.write_text(
                "REMARK 900 GNINA CNN SCORE 0.9274\n"
                "REMARK 900 POSE CONFIDENCE high\n"
                "REMARK 900 POCKET ID "
                "fpocket_02_pocket_2\n"
                "ATOM      1  C   LIG A   1       "
                "0.000   0.000   0.000\n",
                encoding="utf-8",
            )

            hypotheses = collect_pdb_hypotheses(
                root
            )

            self.assertEqual(
                len(hypotheses),
                1,
            )
            self.assertEqual(
                hypotheses[0]["compound"],
                "darunavir",
            )
            self.assertEqual(
                hypotheses[0]["pocket_id"],
                "fpocket_02_pocket_2",
            )
            self.assertAlmostEqual(
                hypotheses[0]["gnina_cnn_score"],
                0.9274,
            )

    def test_write_run_report_combines_target_and_docking_outputs(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)

            target_evidence = {
                "source": {
                    "status": "complete",
                    "result_counts": {
                        "cdd": 1,
                        "interpro": 1,
                        "vogdb": 10,
                    },
                    "tool_statuses": {
                        "cdd": "complete",
                        "interpro": "complete",
                        "vogdb": "complete",
                    },
                    "tool_errors": {},
                },
                "target_interpretation": {
                    "target_name": (
                        "HIV-like retroviral "
                        "aspartyl protease"
                    ),
                    "target_class": (
                        "viral protease"
                    ),
                    "enzyme_class": (
                        "aspartyl protease"
                    ),
                    "viral_family": (
                        "Retroviridae-like / "
                        "retroviral"
                    ),
                    "evidence_confidence": "high",
                    "docking_priority": "high",
                },
                "evidence": {
                    "confidence_reasoning": [
                        (
                            "Specific CDD/Pfam "
                            "retroviral protease "
                            "domain evidence was "
                            "detected."
                        )
                    ],
                    "special_domain_evidence": {
                        "label": (
                            "Retroviral aspartyl "
                            "protease domain"
                        ),
                        "tool": "cdd",
                        "hit_name": "RVP",
                        "accession": "pfam00077",
                        "start": "5",
                        "end": "98",
                        "evalue": "5.7e-36",
                        "score": "117.4",
                    },
                },
                "future_ligand_database_query": {
                    "query_terms": [
                        "HIV protease inhibitor",
                        (
                            "viral aspartyl protease "
                            "inhibitor"
                        ),
                    ]
                },
                "limitations": [
                    (
                        "This target evidence packet "
                        "is generated from "
                        "computational annotation only."
                    )
                ],
            }

            (
                root
                / "target_evidence.json"
            ).write_text(
                json.dumps(target_evidence),
                encoding="utf-8",
            )

            hypothesis_path = (
                root
                / (
                    "01__darunavir__"
                    "fpocket_02_pocket_2__"
                    "hypothesis_01.pdb"
                )
            )

            hypothesis_path.write_text(
                "REMARK 900 GNINA CNN SCORE 0.9274\n"
                "REMARK 900 POSE CONFIDENCE moderate\n"
                "REMARK 900 POCKET ID "
                "fpocket_02_pocket_2\n",
                encoding="utf-8",
            )

            report_path = write_run_report(
                output_dir=root
            )
            report = report_path.read_text(
                encoding="utf-8"
            )

            self.assertIn(
                (
                    "HIV-like retroviral "
                    "aspartyl protease"
                ),
                report,
            )
            self.assertIn(
                "HIV protease inhibitor",
                report,
            )
            self.assertIn(
                "darunavir",
                report,
            )
            self.assertIn(
                "0.9274",
                report,
            )
            self.assertIn(
                "computational target annotation",
                report,
            )
            self.assertIn(
                "| InterPro | complete | 1 |",
                report,
            )

    def test_main_report_preserves_failed_interpro_status(
        self,
    ) -> None:
        error = (
            "RuntimeError: [InterPro] Command failed\n"
            "Exit code: 1\n"
            "Error: bad file format in HMM file "
            "19.0/ncbifam.hmm"
        )

        target_evidence = {
            "schema_version": (
                "target_evidence.v0.2"
            ),
            "source": {
                "job_id": "borna-regression",
                "status": "partial",
                "source_fasta": (
                    "BoDV1_test.fasta"
                ),
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
            },
            "target_interpretation": {
                "target_name": (
                    "viral polymerase"
                ),
                "target_class": (
                    "viral polymerase"
                ),
                "enzyme_class": "polymerase",
                "viral_family": "unknown",
                "evidence_confidence": "medium",
                "docking_priority": "high",
            },
            "evidence": {
                "confidence_reasoning": [
                    (
                        "Usable annotation rows "
                        "retained: CDD=1, "
                        "InterPro=0, VOGDB=6."
                    ),
                    (
                        "CDD completed successfully "
                        "with 1 usable result(s)."
                    ),
                    (
                        "InterPro failed; 0 usable "
                        "result(s) were retained. "
                        "Error summary: Error: bad "
                        "file format in HMM file "
                        "19.0/ncbifam.hmm"
                    ),
                    (
                        "VOGDB completed successfully "
                        "with 6 usable result(s)."
                    ),
                ],
                "special_domain_evidence": None,
            },
            "future_ligand_database_query": {
                "query_terms": [
                    "viral polymerase inhibitor",
                    (
                        "RNA dependent RNA "
                        "polymerase inhibitor"
                    ),
                ],
            },
            "limitations": [
                (
                    "This target evidence packet "
                    "is generated from computational "
                    "annotation only."
                ),
                (
                    "InterPro failed; its zero usable "
                    "results must not be interpreted "
                    "as a successful no-hit result. "
                    "Error summary: Error: bad file "
                    "format in HMM file "
                    "19.0/ncbifam.hmm"
                ),
            ],
        }

        section = "\n".join(
            _render_target_section(
                target_evidence
            )
        )

        self.assertIn(
            (
                "- Overall CPU annotation "
                "status: partial"
            ),
            section,
        )
        self.assertIn(
            "| CDD | complete | 1 |",
            section,
        )
        self.assertIn(
            "| InterPro | failed | 0 |",
            section,
        )
        self.assertIn(
            "| VOGDB | complete | 6 |",
            section,
        )
        self.assertIn(
            (
                "bad file format in HMM file "
                "19.0/ncbifam.hmm"
            ),
            section,
        )
        self.assertIn(
            (
                "Their zero usable-result counts "
                "must not be interpreted as "
                "successful no-hit results."
            ),
            section,
        )
        self.assertIn(
            "### Annotation Limitations",
            section,
        )
        self.assertIn(
            "InterPro failed",
            section,
        )

        self.assertNotIn(
            "| InterPro | complete_no_hits | 0 |",
            section,
        )
        self.assertNotIn(
            (
                "InterPro completed successfully "
                "and returned no hits"
            ),
            section,
        )

    def test_legacy_evidence_reports_unknown_tool_status(
        self,
    ) -> None:
        target_evidence = {
            "source": {
                "status": "partial_or_failed",
                "result_counts": {
                    "cdd": 1,
                    "interpro": 0,
                    "vogdb": 6,
                },
            },
            "target_interpretation": {
                "target_name": (
                    "viral polymerase"
                ),
                "target_class": (
                    "viral polymerase"
                ),
                "enzyme_class": "polymerase",
                "viral_family": "unknown",
                "evidence_confidence": "medium",
                "docking_priority": "high",
            },
            "evidence": {
                "confidence_reasoning": [],
            },
            "future_ligand_database_query": {
                "query_terms": [],
            },
            "limitations": [],
        }

        section = "\n".join(
            _render_target_section(
                target_evidence
            )
        )

        self.assertIn(
            (
                "Per-tool execution statuses "
                "were not available"
            ),
            section,
        )
        self.assertIn(
            (
                "Result counts: "
                "{'cdd': 1, 'interpro': 0, "
                "'vogdb': 6}"
            ),
            section,
        )

        # Legacy counts remain visible, but the report
        # must not invent successful no-hit statuses.
        self.assertNotIn(
            "complete_no_hits",
            section,
        )


if __name__ == "__main__":
    unittest.main()
