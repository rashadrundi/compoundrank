import unittest

from compoundrank.chembl_search import (
    build_chembl_target_terms,
    make_chembl_candidate,
    retrieve_chembl_candidates,
    search_chembl_targets,
    search_chembl_targets_sequence_first,
)


class ChemblSearchTests(unittest.TestCase):
    def test_target_terms_strip_ligand_suffixes_and_deduplicate(self):
        queries = [
            {
                "query": "Example viral protease inhibitor",
                "specificity": 95,
                "retrieval_route": "generic_exact_target_search",
            },
            {
                "query": "Example viral protease ligand",
                "specificity": 90,
                "retrieval_route": "generic_exact_target_search",
            },
            {
                "query": "pfam99999",
                "specificity": 82,
                "retrieval_route": "generic_domain_accession_search",
            },
            {
                "query": "viral protease inhibitor",
                "specificity": 55,
                "retrieval_route": "generic_target_class_search",
            },
        ]

        terms = build_chembl_target_terms(queries)
        texts = [item["term"] for item in terms]

        self.assertEqual(
            texts,
            [
                "Example viral protease",
                "viral protease",
            ],
        )
        self.assertNotIn("pfam99999", texts)

    def test_retrieval_normalizes_and_deduplicates_activities(self):
        queries = [
            {
                "query": "Example viral enzyme inhibitor",
                "specificity": 95,
                "retrieval_route": "generic_exact_target_search",
            }
        ]

        def fake_request(resource, params, timeout_seconds):
            self.assertGreater(timeout_seconds, 0)

            if resource == "target/search":
                return {
                    "targets": [
                        {
                            "target_chembl_id": "CHEMBL_FAKE_TARGET",
                            "pref_name": "Example viral enzyme",
                            "organism": "Example virus",
                            "target_type": "SINGLE PROTEIN",
                        }
                    ]
                }

            if resource == "activity":
                self.assertEqual(
                    params["target_chembl_id"],
                    "CHEMBL_FAKE_TARGET",
                )
                return {
                    "activities": [
                        {
                            "activity_id": 101,
                            "assay_chembl_id": "CHEMBL_FAKE_ASSAY_1",
                            "document_chembl_id": "CHEMBL_FAKE_DOC",
                            "molecule_chembl_id": "CHEMBL_FAKE_A",
                            "molecule_pref_name": "Candidate Alpha",
                            "canonical_smiles": "CCO",
                            "standard_type": "IC50",
                            "standard_relation": "=",
                            "standard_value": "10",
                            "standard_units": "nM",
                            "pchembl_value": "8.0",
                            "assay_type": "B",
                            "potential_duplicate": 0,
                            "data_validity_comment": None,
                        },
                        {
                            "activity_id": 102,
                            "assay_chembl_id": "CHEMBL_FAKE_ASSAY_2",
                            "document_chembl_id": "CHEMBL_FAKE_DOC",
                            "molecule_chembl_id": "CHEMBL_FAKE_A",
                            "molecule_pref_name": "Candidate Alpha",
                            "canonical_smiles": "CCO",
                            "standard_type": "Ki",
                            "standard_relation": "=",
                            "standard_value": "100",
                            "standard_units": "nM",
                            "pchembl_value": "7.0",
                            "assay_type": "B",
                            "potential_duplicate": 0,
                            "data_validity_comment": None,
                        },
                        {
                            "activity_id": 103,
                            "molecule_chembl_id": "CHEMBL_REJECTED",
                            "molecule_pref_name": "Rejected Duplicate",
                            "canonical_smiles": "CCC",
                            "standard_type": "IC50",
                            "standard_relation": "=",
                            "standard_value": "1",
                            "standard_units": "nM",
                            "pchembl_value": "9.0",
                            "potential_duplicate": 1,
                            "data_validity_comment": None,
                        },
                    ]
                }

            raise AssertionError(f"Unexpected resource: {resource}")

        candidates, trace = retrieve_chembl_candidates(
            queries,
            max_candidates=10,
            request_json=fake_request,
        )

        self.assertEqual(len(candidates), 1)

        candidate = candidates[0]
        self.assertEqual(
            candidate["compound_name"],
            "Candidate Alpha",
        )
        self.assertEqual(
            candidate["chembl_molecule_id"],
            "CHEMBL_FAKE_A",
        )
        self.assertEqual(candidate["pchembl_value"], 8.0)
        self.assertFalse(candidate["hardcoded_seed"])
        self.assertEqual(
            candidate["discovery_source"],
            "chembl_activity",
        )
        self.assertEqual(
            len(candidate["supporting_activities"]),
            2,
        )

        self.assertEqual(trace["target_count"], 1)
        self.assertEqual(trace["activity_count"], 3)
        self.assertEqual(
            trace["accepted_activity_count"],
            2,
        )
        self.assertEqual(trace["candidate_count"], 1)

    def test_viral_context_rejects_human_aspartyl_targets(self):
        queries = [
            {
                "query": (
                    "HIV-like retroviral aspartyl "
                    "protease inhibitor"
                ),
                "specificity": 95,
                "retrieval_route": (
                    "generic_exact_target_search"
                ),
            }
        ]

        target_context = {
            "target_name": (
                "HIV-like retroviral aspartyl protease"
            ),
            "target_class": "viral protease",
            "enzyme_class": "aspartyl protease",
            "viral_family": (
                "Retroviridae-like / retroviral"
            ),
            "special_domain_label": (
                "Retroviral aspartyl protease domain"
            ),
        }

        activity_target_ids = []

        def fake_request(resource, params, timeout_seconds):
            if resource == "target/search":
                return {
                    "targets": [
                        {
                            "target_chembl_id": "CHEMBL_HUMAN_BACE",
                            "pref_name": "Beta-secretase 1",
                            "organism": "Homo sapiens",
                            "target_type": "SINGLE PROTEIN",
                        },
                        {
                            "target_chembl_id": "CHEMBL_VIRAL_PROTEASE",
                            "pref_name": (
                                "Human immunodeficiency virus "
                                "type 1 protease"
                            ),
                            "organism": (
                                "Human immunodeficiency virus 1"
                            ),
                            "target_type": "SINGLE PROTEIN",
                        },
                    ]
                }

            if resource == "activity":
                target_id = params["target_chembl_id"]
                activity_target_ids.append(target_id)

                return {
                    "activities": [
                        {
                            "activity_id": 501,
                            "assay_chembl_id": "CHEMBL_ASSAY_501",
                            "document_chembl_id": "CHEMBL_DOC_501",
                            "molecule_chembl_id": "CHEMBL_EXTERNAL_501",
                            "molecule_pref_name": "External Candidate",
                            "canonical_smiles": "CCN",
                            "standard_type": "IC50",
                            "standard_relation": "=",
                            "standard_value": "20",
                            "standard_units": "nM",
                            "pchembl_value": "7.7",
                            "assay_type": "B",
                            "potential_duplicate": 0,
                            "data_validity_comment": None,
                        }
                    ]
                }

            raise AssertionError(
                f"Unexpected resource: {resource}"
            )

        candidates, trace = retrieve_chembl_candidates(
            queries,
            target_context=target_context,
            max_targets=5,
            request_json=fake_request,
        )

        target_ids = {
            target["target_chembl_id"]
            for target in trace["targets"]
        }

        self.assertEqual(
            target_ids,
            {"CHEMBL_VIRAL_PROTEASE"},
        )
        self.assertEqual(
            activity_target_ids,
            ["CHEMBL_VIRAL_PROTEASE"],
        )
        self.assertEqual(len(candidates), 1)
        self.assertFalse(candidates[0]["hardcoded_seed"])

    def test_specific_viral_identity_beats_related_proteases(self):
        queries = [
            {
                "query": (
                    "HIV-like retroviral aspartyl "
                    "protease inhibitor"
                ),
                "specificity": 95,
                "retrieval_route": (
                    "generic_exact_target_search"
                ),
            }
        ]

        context = {
            "target_name": (
                "HIV-like retroviral aspartyl protease"
            ),
            "target_class": "viral protease",
            "enzyme_class": "aspartyl protease",
            "viral_family": (
                "Retroviridae-like / retroviral"
            ),
        }

        def fake_request(resource, params, timeout_seconds):
            self.assertEqual(resource, "target/search")

            return {
                "targets": [
                    {
                        "target_chembl_id": "CHEMBL_GENERIC_HIV",
                        "pref_name": "Protease",
                        "organism": (
                            "Human immunodeficiency virus 1"
                        ),
                        "target_type": "SINGLE PROTEIN",
                    },
                    {
                        "target_chembl_id": "CHEMBL_CANONICAL_HIV",
                        "pref_name": (
                            "Human immunodeficiency virus "
                            "type 1 protease"
                        ),
                        "organism": (
                            "Human immunodeficiency virus 1"
                        ),
                        "target_type": "SINGLE PROTEIN",
                    },
                    {
                        "target_chembl_id": "CHEMBL_NOROVIRUS",
                        "pref_name": "Protease",
                        "organism": "Norovirus",
                        "target_type": "SINGLE PROTEIN",
                    },
                    {
                        "target_chembl_id": "CHEMBL_HTLV",
                        "pref_name": "Protease",
                        "organism": (
                            "Human T-cell leukemia virus type I"
                        ),
                        "target_type": "SINGLE PROTEIN",
                    },
                ]
            }

        targets, terms = search_chembl_targets(
            queries,
            target_context=context,
            request_json=fake_request,
        )

        self.assertGreater(len(terms), 0)
        self.assertEqual(len(targets), 1)
        self.assertEqual(
            targets[0]["target_chembl_id"],
            "CHEMBL_CANONICAL_HIV",
        )

    def test_sequence_resolution_outranks_text_preference(self):
        queries = [
            {
                "query": "Example viral protease inhibitor",
                "specificity": 95,
                "retrieval_route": "generic_exact_target_search",
            }
        ]

        context = {
            "target_name": "Example viral protease",
            "target_class": "viral protease",
            "enzyme_class": "aspartyl protease",
            "viral_family": "Example virus",
        }

        query_sequence = "AAAACCCCGGGG"

        def fake_request(resource, params, timeout_seconds):
            if resource == "target/search":
                return {
                    "targets": [
                        {
                            "target_chembl_id": "CHEMBL_TEXT_FAVORITE",
                            "pref_name": "Example viral protease",
                            "organism": "Example virus",
                            "target_type": "SINGLE PROTEIN",
                        },
                        {
                            "target_chembl_id": "CHEMBL_SEQUENCE_MATCH",
                            "pref_name": "Protease",
                            "organism": "Example virus",
                            "target_type": "SINGLE PROTEIN",
                        },
                    ]
                }

            if resource == "target/CHEMBL_TEXT_FAVORITE":
                return {
                    "target_chembl_id": "CHEMBL_TEXT_FAVORITE",
                    "pref_name": "Example viral protease",
                    "organism": "Example virus",
                    "target_type": "SINGLE PROTEIN",
                    "target_components": [
                        {
                            "component_id": 1,
                            "accession": "FAKE_TEXT",
                            "component_description": "Protease",
                            "sequence": "AAAACCCCGGGA",
                        }
                    ],
                }

            if resource == "target/CHEMBL_SEQUENCE_MATCH":
                return {
                    "target_chembl_id": "CHEMBL_SEQUENCE_MATCH",
                    "pref_name": "Protease",
                    "organism": "Example virus",
                    "target_type": "SINGLE PROTEIN",
                    "target_components": [
                        {
                            "component_id": 2,
                            "accession": "FAKE_EXACT",
                            "component_description": "Protease",
                            "sequence": "AAAACCCCGGGG",
                        }
                    ],
                }

            raise AssertionError(
                f"Unexpected resource: {resource}"
            )

        targets, search_terms = (
            search_chembl_targets_sequence_first(
                queries,
                target_context=context,
                query_sequence=query_sequence,
                request_json=fake_request,
            )
        )

        self.assertGreater(len(search_terms), 0)
        self.assertEqual(len(targets), 1)
        self.assertEqual(
            targets[0]["target_chembl_id"],
            "CHEMBL_SEQUENCE_MATCH",
        )
        self.assertEqual(
            targets[0]["matched_component_accession"],
            "FAKE_EXACT",
        )
        self.assertEqual(
            targets[0]["sequence_identity"],
            1.0,
        )
        self.assertEqual(
            targets[0]["sequence_query_coverage"],
            1.0,
        )
        self.assertEqual(
            targets[0]["target_resolution_route"],
            "sequence_alignment",
        )


    def test_candidate_separates_submitted_and_reference_targets(
        self,
    ):
        activity = {
            "activity_id": 7001,
            "assay_chembl_id": "CHEMBL_ASSAY_7001",
            "document_chembl_id": "CHEMBL_DOC_7001",
            "molecule_chembl_id": "CHEMBL_MOL_7001",
            "molecule_pref_name": "Example inhibitor",
            "canonical_smiles": "CCO",
            "standard_type": "Ki",
            "standard_relation": "=",
            "standard_value": "1.0",
            "standard_units": "nM",
            "pchembl_value": "9.0",
            "assay_type": "B",
            "assay_description": (
                "Biochemical inhibition assay"
            ),
            "bao_label": (
                "single protein format"
            ),
            "bao_endpoint": "BAO_0000190",
            "confidence_score": 9,
            "potential_duplicate": 0,
            "data_validity_comment": None,
        }

        reference_target = {
            "target_chembl_id": (
                "CHEMBL_REFERENCE"
            ),
            "pref_name": (
                "Replicase polyprotein 1ab"
            ),
            "target_type": "SINGLE PROTEIN",
            "organism": "Example coronavirus",
            "tax_id": 1234,
            "target_resolution_route": (
                "sequence_alignment"
            ),
            "target_search_term": (
                "coronavirus main protease"
            ),
            "target_match_score": 215.0,
            "sequence_identity": 1.0,
            "sequence_query_coverage": 1.0,
            "matched_component_accession": (
                "EXAMPLE_MPRO"
            ),
            "matched_component_description": (
                "Main protease"
            ),
        }

        submitted_context = {
            "target_name": (
                "coronavirus 3C-like "
                "main protease"
            ),
            "target_class": (
                "viral protease"
            ),
            "enzyme_class": (
                "cysteine protease"
            ),
            "viral_family": (
                "Coronaviridae-like / "
                "coronavirus"
            ),
            "confidence": "high",
            "special_domain_label": (
                "Coronavirus 3C-like "
                "main protease domain"
            ),
            "special_domain_accession": (
                "cd21666"
            ),
        }

        candidate = make_chembl_candidate(
            activity,
            reference_target,
            target_context=submitted_context,
        )

        self.assertEqual(
            candidate["target_name"],
            (
                "coronavirus 3C-like "
                "main protease"
            ),
        )
        self.assertEqual(
            candidate["enzyme_class"],
            "cysteine protease",
        )
        self.assertEqual(
            candidate[
                "target_evidence_confidence"
            ],
            "high",
        )
        self.assertEqual(
            candidate["submitted_target"][
                "special_domain_accession"
            ],
            "cd21666",
        )
        self.assertEqual(
            candidate["reference_target"][
                "chembl_target_id"
            ],
            "CHEMBL_REFERENCE",
        )
        self.assertEqual(
            candidate["reference_target"][
                "target_name"
            ],
            "Replicase polyprotein 1ab",
        )
        self.assertEqual(
            candidate["reference_target"][
                "sequence_query_coverage"
            ],
            1.0,
        )

        supporting = candidate[
            "supporting_activities"
        ][0]

        self.assertEqual(
            supporting["document_chembl_id"],
            "CHEMBL_DOC_7001",
        )
        self.assertEqual(
            supporting["assay_description"],
            "Biochemical inhibition assay",
        )
        self.assertEqual(
            supporting["confidence_score"],
            9,
        )
        self.assertEqual(
            supporting["evidence_category"],
            "direct_inhibition_constant",
        )


if __name__ == "__main__":
    unittest.main()
