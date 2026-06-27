from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from compoundrank.receptor import (
    _pdb_atom_element,
    _write_pdb2pqr_input,
    prepare_receptor,
)


HEAVY_CARBON = (
    "ATOM      1  CA  ALA A   1"
    "       0.000   0.000   0.000"
    "  1.00  0.00           C"
)

EXPLICIT_HYDROGEN = (
    "ATOM      2  H   ALA A   1"
    "       0.000   1.000   0.000"
    "  1.00  0.00           H"
)

EXPLICIT_DEUTERIUM = (
    "ATOM      3  D   ALA A   1"
    "       0.000   2.000   0.000"
    "  1.00  0.00           D"
)


class ReceptorHydrogenNormalizationTests(
    unittest.TestCase
):
    def test_strips_explicit_hydrogen_and_deuterium(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            source = root / "source.pdb"
            output = root / "heavy.pdb"

            source.write_text(
                "\n".join(
                    [
                        HEAVY_CARBON,
                        EXPLICIT_HYDROGEN,
                        EXPLICIT_DEUTERIUM,
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            removed = (
                _write_pdb2pqr_input(
                    source,
                    output,
                )
            )

            result = output.read_text(
                encoding="utf-8"
            )

            self.assertEqual(
                removed,
                2,
            )

            self.assertIn(
                " CA  ALA",
                result,
            )

            self.assertNotIn(
                " H   ALA",
                result,
            )

            self.assertNotIn(
                " D   ALA",
                result,
            )

    def test_atom_name_fallback_detects_hydrogen(
        self,
    ) -> None:
        hydrogen_without_element = (
            "ATOM      9  HE2 GLU B 119"
            "       1.326  25.249 102.401"
        )

        self.assertEqual(
            _pdb_atom_element(
                hydrogen_without_element
            ),
            "H",
        )

        self.assertEqual(
            _pdb_atom_element(
                HEAVY_CARBON
            ),
            "C",
        )

    def test_prepare_receptor_uses_heavy_atom_input(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            source = (
                root / "hydrogenated.pdb"
            )

            source.write_text(
                "\n".join(
                    [
                        HEAVY_CARBON,
                        EXPLICIT_HYDROGEN,
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            commands: list[
                list[str]
            ] = []

            def fake_run_command(
                command,
                **kwargs,
            ):
                command = [
                    str(value)
                    for value in command
                ]

                commands.append(
                    command
                )

                output_index = (
                    command.index(
                        "--pdb-output"
                    )
                    + 1
                )

                protonated = Path(
                    command[output_index]
                )

                pqr = Path(
                    command[-1]
                )

                protonated.write_text(
                    HEAVY_CARBON
                    + "\nEND\n",
                    encoding="utf-8",
                )

                pqr.write_text(
                    "PQR\n",
                    encoding="utf-8",
                )

                return SimpleNamespace(
                    stdout="",
                    stderr="",
                )

            def fake_meeko(
                executable,
                protonated_pdb,
                output_prefix,
                expected_output,
            ):
                Path(
                    expected_output
                ).write_text(
                    "PDBQT\n",
                    encoding="utf-8",
                )

            with (
                patch(
                    "compoundrank.receptor."
                    "resolve_executable",
                    side_effect=lambda value, label: (
                        value
                    ),
                ),
                patch(
                    "compoundrank.receptor."
                    "run_command",
                    side_effect=(
                        fake_run_command
                    ),
                ),
                patch(
                    "compoundrank.receptor."
                    "_run_meeko_receptor",
                    side_effect=fake_meeko,
                ),
            ):
                prepared = prepare_receptor(
                    source,
                    root / "cache",
                    pdb2pqr_bin="pdb2pqr",
                    meeko_receptor_bin=(
                        "mk_prepare_receptor.py"
                    ),
                )

            self.assertEqual(
                len(commands),
                1,
            )

            pdb2pqr_input = Path(
                commands[0][-2]
            )

            normalized = (
                pdb2pqr_input.read_text(
                    encoding="utf-8"
                )
            )

            copied_source = (
                prepared.prepared_pdbqt
                .parent
                / "receptor_source.pdb"
            )

            preserved = (
                copied_source.read_text(
                    encoding="utf-8"
                )
            )

            self.assertNotIn(
                " H   ALA",
                normalized,
            )

            self.assertIn(
                " H   ALA",
                preserved,
            )

            self.assertTrue(
                prepared.prepared_pdbqt
                .is_file()
            )


if __name__ == "__main__":
    unittest.main()
