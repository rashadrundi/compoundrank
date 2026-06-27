from __future__ import annotations

import inspect
import tempfile
import unittest
from dataclasses import fields
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from compoundrank.gnina import (
    run_gnina_seed,
)
from compoundrank.models import (
    PocketDefinition,
    PoseRecord,
)


class FakeMolecule:
    def __init__(self) -> None:
        self.properties: dict[
            str,
            object,
        ] = {}

    def SetProp(
        self,
        name: str,
        value: str,
    ) -> None:
        self.properties[name] = value

    def SetIntProp(
        self,
        name: str,
        value: int,
    ) -> None:
        self.properties[name] = value

    def SetDoubleProp(
        self,
        name: str,
        value: float,
    ) -> None:
        self.properties[name] = value

    def HasProp(
        self,
        name: str,
    ) -> bool:
        return name in self.properties

    def GetProp(
        self,
        name: str,
    ) -> str:
        return str(
            self.properties[name]
        )


def receptor() -> SimpleNamespace:
    return SimpleNamespace(
        source_pdb=Path(
            "/tmp/snapshot.pdb"
        ),
        prepared_pdbqt=Path(
            "/tmp/snapshot.pdbqt"
        ),
        display_pdb=Path(
            "/tmp/snapshot_display.pdb"
        ),
        cache_key="snapshot",
    )


def ligand() -> SimpleNamespace:
    return SimpleNamespace(
        name="ligand",
        source_sdf=Path(
            "/tmp/ligand.sdf"
        ),
        prepared_pdbqt=Path(
            "/tmp/ligand.pdbqt"
        ),
        cache_key="ligand",
    )


def pocket() -> PocketDefinition:
    return PocketDefinition(
        mode="explicit",
        center_x=0.0,
        center_y=0.0,
        center_z=0.0,
        size_x=20.0,
        size_y=20.0,
        size_z=20.0,
        pocket_id="pocket_01",
    )


class GninaConformerAttributionTests(
    unittest.TestCase
):
    def test_pose_record_has_conformer_fields(
        self,
    ) -> None:
        names = {
            field.name
            for field in fields(
                PoseRecord
            )
        }

        self.assertIn(
            "receptor_conformer_id",
            names,
        )

        self.assertIn(
            "receptor_source_pdb",
            names,
        )

        self.assertIn(
            "receptor_display_pdb",
            names,
        )

    def test_seed_signature_has_conformer_id(
        self,
    ) -> None:
        parameter = (
            inspect.signature(
                run_gnina_seed
            )
            .parameters[
                "receptor_conformer_id"
            ]
        )

        self.assertEqual(
            parameter.default,
            "submitted_receptor",
        )

    def test_rejects_empty_conformer_id(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                ValueError,
                "cannot be empty",
            ):
                run_gnina_seed(
                    receptor(),
                    ligand(),
                    pocket(),
                    1,
                    Path(directory),
                    receptor_conformer_id="",
                )

    def test_rejects_path_separator(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                ValueError,
                "path separators",
            ):
                run_gnina_seed(
                    receptor(),
                    ligand(),
                    pocket(),
                    1,
                    Path(directory),
                    receptor_conformer_id=(
                        "../snapshot"
                    ),
                )

    @patch(
        "compoundrank.gnina."
        "reconstruct_heavy_pose"
    )
    @patch(
        "compoundrank.gnina."
        "choose_pose_to_source_mapping"
    )
    @patch(
        "compoundrank.gnina."
        "parse_meeko_index_pairs"
    )
    @patch(
        "compoundrank.gnina."
        "load_first_sdf"
    )
    @patch(
        "compoundrank.gnina."
        "load_sdf_records"
    )
    @patch(
        "compoundrank.gnina."
        "run_command"
    )
    @patch(
        "compoundrank.gnina."
        "resolve_executable"
    )
    def test_records_conformer_and_isolates_path(
        self,
        resolve_executable,
        run_command,
        load_sdf_records,
        load_first_sdf,
        parse_pairs,
        choose_mapping,
        reconstruct,
    ) -> None:
        raw_pose = FakeMolecule()
        reconstructed = FakeMolecule()

        resolve_executable.return_value = (
            "gnina"
        )

        def fake_run_command(
            command,
            **kwargs,
        ):
            output_path = Path(
                command[
                    command.index(
                        "--out"
                    )
                    + 1
                ]
            )

            output_path.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

            output_path.write_text(
                "placeholder\n",
                encoding="utf-8",
            )

            return SimpleNamespace(
                stdout="",
                stderr="",
            )

        run_command.side_effect = (
            fake_run_command
        )

        load_sdf_records.return_value = [
            raw_pose
        ]

        load_first_sdf.return_value = (
            object()
        )

        parse_pairs.return_value = []
        choose_mapping.return_value = {}
        reconstruct.return_value = (
            reconstructed
        )

        with tempfile.TemporaryDirectory() as directory:
            work_dir = Path(directory)

            records = run_gnina_seed(
                receptor(),
                ligand(),
                pocket(),
                17,
                work_dir,
                receptor_conformer_id=(
                    "snapshot_0001"
                ),
            )

            self.assertEqual(
                len(records),
                1,
            )

            record = records[0]

            self.assertEqual(
                record.receptor_conformer_id,
                "snapshot_0001",
            )

            self.assertEqual(
                record.receptor_source_pdb,
                receptor().source_pdb,
            )

            self.assertEqual(
                record.receptor_display_pdb,
                receptor().display_pdb,
            )

            self.assertEqual(
                reconstructed.properties[
                    "receptor_conformer_id"
                ],
                "snapshot_0001",
            )

            expected = (
                work_dir
                / "snapshot_0001"
                / "ligand"
                / "pocket_01"
                / "seed_17"
                / "poses.sdf"
            )

            self.assertEqual(
                record.source_sdf,
                expected,
            )


if __name__ == "__main__":
    unittest.main()
