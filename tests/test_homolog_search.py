from __future__ import annotations

import unittest

from compoundrank.homolog_search import parse_cpu_response


class HomologSearchTests(unittest.TestCase):
    def test_parse_cpu_response_counts_rows(self) -> None:
        response = {
            "job_id": "abc",
            "status": "complete",
            "results": {
                "cdd": [{"id": 1}, {"id": 2}],
                "interpro": [{"id": 3}],
                "vogdb": [],
            },
            "files": {
                "report": "example.json",
            },
        }

        parsed = parse_cpu_response(response)

        self.assertEqual(parsed["job_id"], "abc")
        self.assertEqual(parsed["status"], "complete")
        self.assertEqual(parsed["result_counts"]["cdd"], 2)
        self.assertEqual(parsed["result_counts"]["interpro"], 1)
        self.assertEqual(parsed["result_counts"]["vogdb"], 0)
        self.assertEqual(parsed["files"]["report"], "example.json")


if __name__ == "__main__":
    unittest.main()
