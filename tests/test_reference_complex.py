import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.reference_complex import extract_reference_complex, summarize_ligands


def pdb_line(record, serial, atom, resname, chain, resseq, x, y, z, element="C", altloc=" "):
    return (
        f"{record:<6}{serial:>5} "
        f"{atom:<4}{altloc}{resname:>3} {chain:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"  1.00 20.00          {element:>2}"
    )


class ReferenceComplexTests(unittest.TestCase):
    def test_summarize_and_extract_reference_ligand(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            complex_pdb = tmp / "complex.pdb"
            complex_pdb.write_text(
                "\n".join(
                    [
                        pdb_line("ATOM", 1, "N", "PRO", "A", 1, 0.0, 0.0, 0.0, "N"),
                        pdb_line("ATOM", 2, "CA", "PRO", "A", 1, 1.0, 0.0, 0.0, "C"),
                        pdb_line("HETATM", 3, "C1", "017", "A", 501, 10.0, 11.0, 12.0, "C"),
                        pdb_line("HETATM", 4, "C2", "017", "A", 501, 12.0, 13.0, 14.0, "C"),
                        pdb_line("HETATM", 5, "O", "HOH", "A", 900, 99.0, 99.0, 99.0, "O"),
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            summaries = summarize_ligands(complex_pdb)
            self.assertEqual(len(summaries), 1)
            self.assertEqual(summaries[0].resname, "017")
            self.assertEqual(summaries[0].atom_count, 2)
            self.assertAlmostEqual(summaries[0].center_x, 11.0)

            output_dir = tmp / "extract"
            outputs = extract_reference_complex(
                complex_pdb=complex_pdb,
                output_dir=output_dir,
                ligand_resname="017",
                ligand_chain="A",
                ligand_resseq="501",
                padding=8.0,
                min_size=18.0,
            )

            self.assertTrue(outputs["reference_receptor"].exists())
            self.assertTrue(outputs["reference_ligand"].exists())
            self.assertTrue(outputs["reference_box"].exists())

            ligand_text = outputs["reference_ligand"].read_text(encoding="utf-8")
            self.assertIn("017", ligand_text)
            self.assertNotIn("HOH", ligand_text)

            box = json.loads(outputs["reference_box"].read_text(encoding="utf-8"))
            self.assertEqual(box["box_mode"], "reference_ligand")
            self.assertAlmostEqual(box["center_x"], 11.0)
            self.assertAlmostEqual(box["center_y"], 12.0)
            self.assertAlmostEqual(box["center_z"], 13.0)


if __name__ == "__main__":
    unittest.main()
