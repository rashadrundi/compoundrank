from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from compoundrank.pocket import parse_fpocket_info


class FpocketTopNTests(unittest.TestCase):
    def test_fpocket_parser_reads_multiple_pockets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            info = Path(temporary) / "fpocket_info.txt"
            info.write_text(
                "Pocket 1:\n\tScore: 0.15\nPocket 2:\n\tScore: 0.42\n",
                encoding="utf-8",
            )
            pockets = parse_fpocket_info(info)

        self.assertEqual(len(pockets), 2)
        self.assertEqual(pockets[1]["number"], 2)
        self.assertEqual(pockets[1]["metrics"]["Score"], 0.42)


if __name__ == "__main__":
    unittest.main()
