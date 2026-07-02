import unittest

from compoundrank.target_evidence import TARGET_RULES


class TargetEvidenceHelicaseSpecificityTests(unittest.TestCase):
    def test_generic_ntpase_atpase_are_weak_not_primary_helicase_terms(self):
        helicase_rules = [
            rule
            for rule in TARGET_RULES
            if rule.get("target_class") == "viral helicase"
        ]

        self.assertTrue(helicase_rules)

        rule = helicase_rules[0]
        keywords = {
            str(term).lower()
            for term in rule.get("keywords", [])
        }
        weak_keywords = {
            str(term).lower()
            for term in rule.get("weak_keywords", [])
        }

        self.assertIn("helicase", keywords)
        self.assertNotIn("ntpase", keywords)
        self.assertNotIn("atpase", keywords)
        self.assertNotIn("nucleoside triphosphatase", keywords)

        self.assertIn("ntpase", weak_keywords)
        self.assertIn("atpase", weak_keywords)
        self.assertIn("nucleoside triphosphatase", weak_keywords)


if __name__ == "__main__":
    unittest.main()
