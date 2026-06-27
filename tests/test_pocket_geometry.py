from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.models import PocketDefinition
from compoundrank.pocket import (
    _box_from_pocket_file,
    _coordinate_file_for_pocket,
    write_pocket_definitions,
)


def atom_line(
    serial: int,
    x: float,
    y: float,
    z: float,
) -> str:
    return (
        f"ATOM  {serial:5d}  C   UNK A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        "  1.00  0.00           C\n"
    )


class PocketGeometryTests(unittest.TestCase):
    def test_vertex_file_is_preferred_over_atom_file(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            out_dir = Path(directory)
            pockets = out_dir / "pockets"
            pockets.mkdir()

            atom_path = pockets / "pocket1_atm.pdb"
            vertex_path = (
                pockets / "pocket1_vert.pqr"
            )

            atom_path.write_text(
                atom_line(1, 100.0, 100.0, 100.0),
                encoding="utf-8",
            )
            vertex_path.write_text(
                atom_line(1, 0.0, 0.0, 0.0),
                encoding="utf-8",
            )

            selected = _coordinate_file_for_pocket(
                out_dir,
                1,
            )

            self.assertEqual(
                selected,
                vertex_path,
            )

    def test_box_uses_twenty_angstrom_minimum(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pocket1_vert.pqr"

            path.write_text(
                atom_line(1, 0.0, 0.0, 0.0)
                + atom_line(2, 15.0, 2.0, 3.0),
                encoding="utf-8",
            )

            center, size = _box_from_pocket_file(
                path,
                padding=4.0,
            )

            self.assertAlmostEqual(center[0], 7.5)
            self.assertAlmostEqual(center[1], 1.0)
            self.assertAlmostEqual(center[2], 1.5)

            self.assertAlmostEqual(size[0], 23.0)
            self.assertAlmostEqual(size[1], 20.0)
            self.assertAlmostEqual(size[2], 20.0)

    def test_writes_machine_readable_definitions(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = (
                Path(directory)
                / "pocket_definitions.json"
            )

            pocket = PocketDefinition(
                mode="explicit",
                center_x=1.0,
                center_y=2.0,
                center_z=3.0,
                size_x=20.0,
                size_y=20.0,
                size_z=20.0,
                source="fpocket rank 1",
                pocket_id="fpocket_01_pocket_3",
                pocket_rank=1,
                fpocket_score=0.5,
            )

            write_pocket_definitions(
                output,
                [pocket],
            )

            payload = json.loads(
                output.read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                payload["pocket_count"],
                1,
            )
            self.assertFalse(
                payload[
                    "reference_ligand_used_for_selection"
                ]
            )
            self.assertEqual(
                payload["pockets"][0]["pocket_id"],
                "fpocket_01_pocket_3",
            )


if __name__ == "__main__":
    unittest.main()
