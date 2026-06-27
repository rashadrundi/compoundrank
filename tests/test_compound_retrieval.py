import json
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from compoundrank.compound_retrieval import (
    build_ligand_candidates,
    fetch_candidate_structures,
    run_compound_retrieval,
)


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

    def test_generic_strict_does_not_return_local_seed_compounds(self):
        target_evidence = {
            "target_name": "HIV-like retroviral aspartyl protease",
            "target_class": "viral protease",
            "enzyme_class": "aspartyl protease",
            "evidence_confidence": "high",
            "supporting_hits": [
                {
                    "source": "CDD",
                    "name": "RVP",
                    "accession": "pfam00077",
                }
            ],
        }

        candidates = build_ligand_candidates(
            target_evidence,
            max_candidates=20,
            retrieval_mode="generic-strict",
        )

        self.assertEqual(candidates, [])

    def test_generic_strict_writes_clean_query_plan_and_metadata(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            target_path = root / "target_evidence.json"
            out_dir = root / "strict"

            target_path.write_text(
                json.dumps(
                    {
                        "target_name": "HIV-like retroviral aspartyl protease",
                        "target_class": "viral protease",
                        "enzyme_class": "aspartyl protease",
                        "evidence_confidence": "high",
                        "supporting_hits": [
                            {
                                "source": "CDD",
                                "name": "RVP",
                                "accession": "pfam00077",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            empty_trace = {
                "backend": "ChEMBL",
                "enabled": True,
                "search_terms": [],
                "target_count": 0,
                "targets": [],
                "activity_count": 0,
                "accepted_activity_count": 0,
                "candidate_count": 0,
            }

            with patch(
                "compoundrank.compound_retrieval."
                "retrieve_chembl_candidates",
                return_value=([], empty_trace),
            ):
                outputs = run_compound_retrieval(
                    target_evidence_path=target_path,
                    output_dir=out_dir,
                    fetch_structures=False,
                    retrieval_mode="generic-strict",
                )

            metadata = json.loads(
                outputs["retrieval_metadata"].read_text(encoding="utf-8")
            )
            query_plan = json.loads(
                outputs["generic_search_queries"].read_text(encoding="utf-8")
            )
            candidate_payload = json.loads(
                outputs["ligand_candidates"].read_text(encoding="utf-8")
            )

            self.assertEqual(
                metadata["retrieval_mode"],
                "generic-strict",
            )
            self.assertFalse(
                metadata["local_rule_registry_enabled"]
            )
            self.assertEqual(
                metadata["hardcoded_candidates_used"],
                0,
            )
            self.assertTrue(
                metadata["strict_provenance_passed"]
            )
            self.assertEqual(
                candidate_payload["candidates"],
                [],
            )
            self.assertGreater(
                len(query_plan["queries"]),
                0,
            )

            serialized_queries = json.dumps(query_plan).lower()
            self.assertNotIn("darunavir", serialized_queries)
            self.assertNotIn("saquinavir", serialized_queries)

    def test_generic_strict_accepts_external_chembl_candidates(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            target_path = root / "target_evidence.json"
            out_dir = root / "strict_external"

            target_path.write_text(
                json.dumps(
                    {
                        "target_name": "Example viral enzyme",
                        "target_class": "viral enzyme",
                        "enzyme_class": "hydrolase",
                        "evidence_confidence": "high",
                    }
                ),
                encoding="utf-8",
            )

            external_candidate = {
                "compound_name": "Candidate Alpha",
                "retrieval_rank": 1,
                "design_status": "database_observed_ligand",
                "discovery_source": "chembl_activity",
                "hardcoded_seed": False,
                "retrieval_route": (
                    "generic_chembl_target_activity_search"
                ),
                "retrieval_rule_id": None,
                "retrieval_reason": "Mock measured ChEMBL activity.",
                "target_family_basis": "Example viral enzyme",
                "target_name": "Example viral enzyme",
                "target_class": "SINGLE PROTEIN",
                "enzyme_class": None,
                "viral_family": "Example virus",
                "special_domain_label": None,
                "special_domain_accession": None,
                "evidence_level": "strong",
                "source_databases": ["ChEMBL"],
                "chembl_target_id": "CHEMBL_FAKE_TARGET",
                "chembl_molecule_id": "CHEMBL_FAKE_A",
                "chembl_activity_id": 101,
                "pchembl_value": 8.0,
                "pubchem_cid": None,
                "structure_fetch_status": "not_attempted",
                "local_sdf_path": None,
                "selected_for_docking": True,
            }

            trace = {
                "backend": "ChEMBL",
                "enabled": True,
                "search_terms": [],
                "target_count": 1,
                "targets": [],
                "activity_count": 1,
                "accepted_activity_count": 1,
                "candidate_count": 1,
            }

            with patch(
                "compoundrank.compound_retrieval."
                "retrieve_chembl_candidates",
                return_value=([external_candidate], trace),
            ):
                outputs = run_compound_retrieval(
                    target_evidence_path=target_path,
                    output_dir=out_dir,
                    fetch_structures=False,
                    retrieval_mode="generic-strict",
                )

            payload = json.loads(
                outputs["ligand_candidates"].read_text(
                    encoding="utf-8"
                )
            )
            metadata = json.loads(
                outputs["retrieval_metadata"].read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(len(payload["candidates"]), 1)
            candidate = payload["candidates"][0]

            self.assertEqual(
                candidate["compound_name"],
                "Candidate Alpha",
            )
            self.assertEqual(
                candidate["discovery_source"],
                "chembl_activity",
            )
            self.assertFalse(candidate["hardcoded_seed"])

            self.assertEqual(
                metadata["external_database_backends_enabled"],
                ["ChEMBL"],
            )
            self.assertEqual(
                metadata["chembl_candidate_count"],
                1,
            )
            self.assertEqual(
                metadata["hardcoded_candidates_used"],
                0,
            )
            self.assertTrue(
                metadata["strict_provenance_passed"]
            )

    def test_chembl_structure_fetch_uses_molecule_id(self):
        candidate = {
            "compound_name": "CHEMBL_FAKE_A",
            "chembl_molecule_id": "CHEMBL_FAKE_A",
            "discovery_source": "chembl_activity",
            "hardcoded_seed": False,
            "source_databases": ["ChEMBL"],
            "selected_for_docking": True,
            "structure_fetch_status": "not_attempted",
            "local_sdf_path": None,
        }

        molecule_payload = {
            "molecule_chembl_id": "CHEMBL_FAKE_A",
            "pref_name": "Candidate Alpha",
            "molecule_type": "Small molecule",
            "max_phase": 0,
            "molecule_structures": {
                "canonical_smiles": "CCO",
                "standard_inchi": (
                    "InChI=1S/C2H6O/"
                    "c1-2-3/h3H,2H2,1H3"
                ),
                "standard_inchi_key": (
                    "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
                ),
            },
        }

        sdf_data = (
            b"Candidate Alpha\n"
            b"  CompoundRank\n"
            b"\n"
            + (b" " * 140)
            + b"\nM  END\n$$$$\n"
        )

        download_calls = []

        def fake_download(
            url,
            *,
            accept,
            timeout_seconds,
        ):
            download_calls.append(
                {
                    "url": url,
                    "accept": accept,
                    "timeout_seconds": timeout_seconds,
                }
            )

            if url.endswith(".json"):
                return json.dumps(
                    molecule_payload
                ).encode("utf-8")

            if url.endswith(".sdf"):
                return sdf_data

            raise AssertionError(
                f"Unexpected URL: {url}"
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            with patch(
                "compoundrank.compound_retrieval."
                "_download_url_bytes",
                side_effect=fake_download,
            ):
                with patch(
                    "compoundrank.compound_retrieval."
                    "fetch_candidate_structures_from_pubchem"
                ) as pubchem_fetch:
                    result = fetch_candidate_structures(
                        [candidate],
                        output_dir=Path(tmpdir),
                    )

            hydrated = result[0]

            pubchem_fetch.assert_not_called()

            self.assertEqual(len(download_calls), 2)
            self.assertTrue(
                any(
                    call["url"].endswith(
                        "/CHEMBL_FAKE_A.json"
                    )
                    for call in download_calls
                )
            )
            self.assertTrue(
                any(
                    call["url"].endswith(
                        "/CHEMBL_FAKE_A.sdf"
                    )
                    for call in download_calls
                )
            )

            self.assertEqual(
                hydrated["compound_name"],
                "Candidate Alpha",
            )
            self.assertEqual(
                hydrated["smiles"],
                "CCO",
            )
            self.assertEqual(
                hydrated["inchi_key"],
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            )
            self.assertEqual(
                hydrated["structure_fetch_status"],
                "fetched",
            )
            self.assertEqual(
                hydrated["structure_source"],
                "ChEMBL molecule SDF endpoint",
            )
            self.assertTrue(
                Path(
                    hydrated["local_sdf_path"]
                ).exists()
            )


if __name__ == "__main__":
    unittest.main()
