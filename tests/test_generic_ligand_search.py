import unittest

from compoundrank.generic_ligand_search import generate_generic_queries


class GenericLigandSearchTests(unittest.TestCase):
    def test_generates_queries_without_target_specific_compound_names(self):
        target = {
            "target_name": "Influenza A virus neuraminidase",
            "target_class": "viral neuraminidase",
            "enzyme_class": "glycosidase",
            "special_domain_label": "Influenza neuraminidase",
            "special_domain_accession": "PF00677",
            "query_terms": [],
        }

        queries = generate_generic_queries(target)
        texts = [item["query"] for item in queries]

        self.assertIn(
            "Influenza A virus neuraminidase inhibitor",
            texts,
        )
        self.assertIn(
            "Influenza neuraminidase inhibitor",
            texts,
        )
        self.assertIn("glycosidase inhibitor", texts)
        self.assertIn("PF00677", texts)

        combined = " ".join(texts).lower()
        self.assertNotIn("oseltamivir", combined)
        self.assertNotIn("zanamivir", combined)
        self.assertNotIn("peramivir", combined)
        self.assertNotIn("laninamivir", combined)

    def test_evidence_queries_are_prioritized_and_deduplicated(self):
        target = {
            "target_name": "Influenza neuraminidase",
            "target_class": "viral neuraminidase",
            "enzyme_class": "glycosidase",
            "special_domain_label": "Neuraminidase",
            "special_domain_accession": "",
            "query_terms": [
                "neuraminidase inhibitor",
                "Neuraminidase inhibitor",
                "influenza antiviral ligand",
            ],
        }

        queries = generate_generic_queries(target)
        texts = [item["query"].casefold() for item in queries]

        self.assertEqual(texts.count("neuraminidase inhibitor"), 1)
        self.assertEqual(
            queries[0]["retrieval_route"],
            "generic_evidence_query",
        )
        self.assertEqual(queries[0]["specificity"], 100)

    def test_unknown_context_does_not_generate_garbage_queries(self):
        target = {
            "target_name": "unknown",
            "target_class": "unknown",
            "enzyme_class": "unknown",
            "special_domain_label": "",
            "special_domain_accession": "",
            "query_terms": [],
        }

        self.assertEqual(generate_generic_queries(target), [])

    def test_max_queries_is_respected(self):
        target = {
            "target_name": "Example viral enzyme",
            "target_class": "viral enzyme",
            "enzyme_class": "hydrolase",
            "special_domain_label": "Example catalytic domain",
            "special_domain_accession": "PF99999",
            "query_terms": ["example inhibitor"],
        }

        queries = generate_generic_queries(target, max_queries=3)
        self.assertEqual(len(queries), 3)


if __name__ == "__main__":
    unittest.main()
