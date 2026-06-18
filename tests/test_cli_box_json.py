import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.cli import load_box_json


class CliBoxJsonTests(unittest.TestCase):
    def test_load_box_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "reference_box.json"
            path.write_text(
                json.dumps(
                    {
                        "center_x": 1,
                        "center_y": 2,
                        "center_z": 3,
                        "size_x": 10,
                        "size_y": 11,
                        "size_z": 12,
                    }
                ),
                encoding="utf-8",
            )

            box = load_box_json(path)

            self.assertEqual(box["center_x"], 1.0)
            self.assertEqual(box["center_y"], 2.0)
            self.assertEqual(box["center_z"], 3.0)
            self.assertEqual(box["size_x"], 10.0)
            self.assertEqual(box["size_y"], 11.0)
            self.assertEqual(box["size_z"], 12.0)

    def test_load_box_json_requires_all_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "bad_box.json"
            path.write_text(json.dumps({"center_x": 1}), encoding="utf-8")

            with self.assertRaises(ValueError):
                load_box_json(path)


if __name__ == "__main__":
    unittest.main()
