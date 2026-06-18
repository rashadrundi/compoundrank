import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.compound_retrieval import build_ligand_candidates, run_compound_retrieval


class CompoundRetrievalTests(unittest.TestCase):
    def test_retroviral_aspartyl_protease_rule_returns_known_inhibitors(self):
        target_evidence = {
            "target_name": "HIV-like retroviral aspartyl protease",
            "target_class": "viral protease",
            "enzyme_class": "aspartyl protease",
            "evidence_confidence": "high",
            "special_domain_evidence": {
                "name": "RVP",
                "accession": "pfam00077",
                "description": "Retroviral aspartyl protease domain",
            },
            "future_ligand_database_query_terms": [
                "HIV protease inhibitor",
                "viral aspartyl protease inhibitor",
            ],
            "supporting_hits": [
                {"source": "CDD", "name": "RVP", "accession": "pfam00077"},
                {"source": "InterPro", "name": "Retropepsins"},
            ],
        }

        candidates = build_ligand_candidates(target_evidence, max_candidates=20)
        names = {candidate["compound_name"].lower() for candidate in candidates}

        self.assertIn("darunavir", names)
        self.assertIn("saquinavir", names)

        darunavir = next(candidate for candidate in candidates if candidate["compound_name"].lower() == "darunavir")
        self.assertEqual(darunavir["design_status"], "known_inhibitor")
        self.assertEqual(darunavir["evidence_level"], "strong")
        self.assertEqual(darunavir["retrieval_rule_id"], "retroviral_aspartyl_protease")
        self.assertIn("RVP", darunavir["domain_basis"])
        self.assertIn("pfam00077", darunavir["domain_basis"])

    def test_run_compound_retrieval_writes_outputs_without_network(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            target_path = root / "target_evidence.json"
            out_dir = root / "out"

            target_path.write_text(
                json.dumps(
                    {
                        "target_name": "HIV-like retroviral aspartyl protease",
                        "target_class": "viral protease",
                        "enzyme_class": "aspartyl protease",
                        "evidence_confidence": "high",
                        "supporting_hits": [
                            {"source": "CDD", "name": "RVP", "accession": "pfam00077"},
                            {"source": "InterPro", "name": "Retropepsins"},
                        ],
                    }
                ),
                encoding="utf-8",
            )

            outputs = run_compound_retrieval(
                target_evidence_path=target_path,
                output_dir=out_dir,
                fetch_structures=False,
            )

            for path in outputs.values():
                self.assertTrue(path.exists(), path)

            payload = json.loads(outputs["ligand_candidates"].read_text(encoding="utf-8"))
            names = {candidate["compound_name"].lower() for candidate in payload["candidates"]}
            self.assertIn("darunavir", names)
            self.assertIn("saquinavir", names)

            report = outputs["ligand_search_report"].read_text(encoding="utf-8")
            self.assertIn("retroviral_aspartyl_protease", report)
            self.assertIn("darunavir", report)

    def test_docking_manifest_uses_main_pipeline_schema(self):
        import tempfile
        from pathlib import Path

        from compoundrank.compound_retrieval import write_docking_manifest

        candidate = {
            "compound_name": "darunavir",
            "local_sdf_path": "/tmp/darunavir.sdf",
            "selected_for_docking": True,
            "retrieval_reason": "test reason",
            "evidence_level": "strong",
            "design_status": "known_inhibitor",
            "pubchem_cid": "213039",
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "docking_manifest.csv"
            write_docking_manifest(output, [candidate])

            text = output.read_text(encoding="utf-8")
            first_line = text.splitlines()[0]
            self.assertEqual(
                first_line,
                "name,source_type,value,retrieval_reason,evidence_level,design_status,pubchem_cid",
            )
            self.assertIn("darunavir,file,/tmp/darunavir.sdf", text)


if __name__ == "__main__":
    unittest.main()
