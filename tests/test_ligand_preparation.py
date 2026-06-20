import csv
import json
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from compoundrank.ligand_preparation import (
    prepare_ligand_manifest,
)


class LigandPreparationTests(unittest.TestCase):
    def test_prepares_3d_ligand_and_removes_counterion(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)

            raw_sdf = root / "raw_ligand.sdf"
            manifest = root / "docking_manifest.csv"
            output_dir = root / "prepared"

            # Flexible organic ligand plus sodium counterion.
            molecule = Chem.MolFromSmiles(
                "CCOCCN.[Na+]"
            )
            self.assertIsNotNone(molecule)

            AllChem.Compute2DCoords(molecule)

            writer = Chem.SDWriter(str(raw_sdf))
            writer.write(molecule)
            writer.close()

            with manifest.open(
                "w",
                encoding="utf-8",
                newline="",
            ) as handle:
                csv_writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "name",
                        "source_type",
                        "value",
                    ],
                )
                csv_writer.writeheader()
                csv_writer.writerow(
                    {
                        "name": "Example Ligand",
                        "source_type": "file",
                        "value": str(raw_sdf),
                    }
                )

            outputs = prepare_ligand_manifest(
                input_manifest=manifest,
                output_dir=output_dir,
                random_seed=12345,
            )

            report = json.loads(
                outputs[
                    "ligand_preparation_report"
                ].read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                report["summary"]["prepared_count"],
                1,
            )
            self.assertEqual(
                report["summary"]["failed_count"],
                0,
            )

            entry = report["ligands"][0]

            self.assertEqual(
                entry["status"],
                "prepared",
            )
            self.assertEqual(
                entry["input_fragment_count"],
                2,
            )
            self.assertGreater(
                entry["explicit_hydrogens"],
                0,
            )
            self.assertTrue(
                entry["marked_as_3d"]
            )
            self.assertGreater(
                entry["z_span"],
                0.01,
            )
            self.assertFalse(
                entry["pH_aware_protonation"]
            )

            prepared_path = Path(
                entry["output_path"]
            )
            self.assertTrue(
                prepared_path.exists()
            )

            supplier = Chem.SDMolSupplier(
                str(prepared_path),
                removeHs=False,
                sanitize=True,
            )
            prepared = next(
                molecule
                for molecule in supplier
                if molecule is not None
            )

            self.assertEqual(
                len(
                    Chem.GetMolFrags(
                        prepared,
                        asMols=False,
                    )
                ),
                1,
            )
            self.assertEqual(
                prepared.GetNumConformers(),
                1,
            )
            self.assertTrue(
                prepared.GetConformer().Is3D()
            )

            with outputs[
                "prepared_docking_manifest"
            ].open(
                "r",
                encoding="utf-8",
                newline="",
            ) as handle:
                rows = list(
                    csv.DictReader(handle)
                )

            self.assertEqual(len(rows), 1)
            self.assertEqual(
                rows[0]["source_type"],
                "file",
            )
            self.assertEqual(
                Path(rows[0]["value"]),
                prepared_path,
            )


if __name__ == "__main__":
    unittest.main()
