from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.run_report import collect_pdb_hypotheses, write_run_report


class RunReportTests(unittest.TestCase):
    def test_collect_pdb_hypotheses_from_filename_and_remarks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            pdb = root / "01__darunavir__fpocket_02_pocket_2__hypothesis_01.pdb"
            pdb.write_text(
                "REMARK 900 GNINA CNN SCORE 0.9274\n"
                "REMARK 900 POSE CONFIDENCE high\n"
                "REMARK 900 POCKET ID fpocket_02_pocket_2\n"
                "ATOM      1  C   LIG A   1       0.000   0.000   0.000\n",
                encoding="utf-8",
            )

            hypotheses = collect_pdb_hypotheses(root)

            self.assertEqual(len(hypotheses), 1)
            self.assertEqual(hypotheses[0]["compound"], "darunavir")
            self.assertEqual(hypotheses[0]["pocket_id"], "fpocket_02_pocket_2")
            self.assertAlmostEqual(hypotheses[0]["gnina_cnn_score"], 0.9274)

    def test_write_run_report_combines_target_and_docking_outputs(self) -> None:
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
                },
                "target_interpretation": {
                    "target_name": "HIV-like retroviral aspartyl protease",
                    "target_class": "viral protease",
                    "enzyme_class": "aspartyl protease",
                    "viral_family": "Retroviridae-like / retroviral",
                    "evidence_confidence": "high",
                    "docking_priority": "high",
                },
                "evidence": {
                    "confidence_reasoning": [
                        "Specific CDD/Pfam retroviral protease domain evidence was detected."
                    ],
                    "special_domain_evidence": {
                        "label": "Retroviral aspartyl protease domain",
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
                        "viral aspartyl protease inhibitor",
                    ]
                },
            }

            (root / "target_evidence.json").write_text(
                json.dumps(target_evidence),
                encoding="utf-8",
            )

            (root / "01__darunavir__fpocket_02_pocket_2__hypothesis_01.pdb").write_text(
                "REMARK 900 GNINA CNN SCORE 0.9274\n"
                "REMARK 900 POSE CONFIDENCE moderate\n"
                "REMARK 900 POCKET ID fpocket_02_pocket_2\n",
                encoding="utf-8",
            )

            report_path = write_run_report(output_dir=root)
            report = report_path.read_text(encoding="utf-8")

            self.assertIn("HIV-like retroviral aspartyl protease", report)
            self.assertIn("HIV protease inhibitor", report)
            self.assertIn("darunavir", report)
            self.assertIn("0.9274", report)
            self.assertIn("computational target annotation", report)


if __name__ == "__main__":
    unittest.main()
