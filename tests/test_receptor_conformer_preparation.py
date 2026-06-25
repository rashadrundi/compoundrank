from __future__ import annotations

import inspect
import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.cli import build_parser
from compoundrank.models import (
    PreparedReceptor,
)
from compoundrank.pipeline import (
    run_pipeline,
)
from compoundrank.receptor_conformer_preparation import (
    prepare_receptor_conformers,
    write_receptor_conformer_preparation,
)


def _write(
    path: Path,
    text: str = "content\n",
) -> Path:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        text,
        encoding="utf-8",
    )

    return path


class ReceptorConformerPreparationTests(
    unittest.TestCase
):
    def test_prepares_submitted_and_snapshots(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            submitted = _write(
                root / "submitted.pdb"
            )

            snapshot_1 = _write(
                root / "snapshot_0001.pdb"
            )

            snapshot_2 = _write(
                root / "snapshot_0002.pdb"
            )

            prepared_calls: list[
                Path
            ] = []

            def fake_prepare(
                source_pdb,
                cache_root,
                **kwargs,
            ):
                source = Path(
                    source_pdb
                )

                prepared_calls.append(
                    source
                )

                token = source.stem

                return PreparedReceptor(
                    source_pdb=source,
                    prepared_pdbqt=_write(
                        root
                        / token
                        / "prepared.pdbqt"
                    ),
                    display_pdb=_write(
                        root
                        / token
                        / "display.pdb"
                    ),
                    cache_key=token,
                )

            aligned = {
                "snapshot_count": 2,
                "source_manifest": (
                    str(
                        root / "aligned.json"
                    )
                ),
                "snapshots": [
                    {
                        "conformer_id": (
                            "snapshot_0001"
                        ),
                        "aligned_path": str(
                            snapshot_1
                        ),
                    },
                    {
                        "conformer_id": (
                            "snapshot_0002"
                        ),
                        "aligned_path": str(
                            snapshot_2
                        ),
                    },
                ],
            }

            conformers = (
                prepare_receptor_conformers(
                    submitted_receptor_pdb=(
                        submitted
                    ),
                    aligned_ensemble=aligned,
                    cache_root=(
                        root / "cache"
                    ),
                    ph=7.4,
                    pdb2pqr_bin="pdb2pqr",
                    meeko_receptor_bin=(
                        "mk_prepare_receptor.py"
                    ),
                    prepare_receptor_fn=(
                        fake_prepare
                    ),
                )
            )

            self.assertEqual(
                [
                    conformer_id
                    for (
                        conformer_id,
                        _,
                    ) in conformers
                ],
                [
                    "submitted_receptor",
                    "snapshot_0001",
                    "snapshot_0002",
                ],
            )

            self.assertEqual(
                prepared_calls,
                [
                    submitted.resolve(),
                    snapshot_1.resolve(),
                    snapshot_2.resolve(),
                ],
            )

    def test_writes_preparation_audit(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            def prepared(
                token: str,
            ) -> PreparedReceptor:
                return PreparedReceptor(
                    source_pdb=_write(
                        root
                        / token
                        / "source.pdb"
                    ),
                    prepared_pdbqt=_write(
                        root
                        / token
                        / "prepared.pdbqt"
                    ),
                    display_pdb=_write(
                        root
                        / token
                        / "display.pdb"
                    ),
                    cache_key=token,
                )

            output = (
                write_receptor_conformer_preparation(
                    root
                    / "preparation.json",
                    [
                        (
                            "submitted_receptor",
                            prepared("submitted"),
                        ),
                        (
                            "snapshot_0001",
                            prepared("snapshot"),
                        ),
                    ],
                    aligned_ensemble={
                        "snapshot_count": 1,
                        "source_manifest": (
                            "/tmp/aligned.json"
                        ),
                    },
                )
            )

            payload = json.loads(
                output.read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                payload[
                    "prepared_conformer_count"
                ],
                2,
            )

            self.assertEqual(
                payload[
                    "docking_behavior"
                ],
                "submitted_receptor_only",
            )

            self.assertEqual(
                [
                    row[
                        "conformer_id"
                    ]
                    for row in payload[
                        "conformers"
                    ]
                ],
                [
                    "submitted_receptor",
                    "snapshot_0001",
                ],
            )

    def test_cli_accepts_aligned_manifest(
        self,
    ) -> None:
        arguments = (
            build_parser().parse_args(
                [
                    "--receptor",
                    "/tmp/receptor.pdb",
                    (
                        "--aligned-receptor-"
                        "ensemble-json"
                    ),
                    "/tmp/aligned.json",
                    "--data-root",
                    "/tmp/data",
                ]
            )
        )

        self.assertEqual(
            arguments
            .aligned_receptor_ensemble_json,
            "/tmp/aligned.json",
        )

    def test_pipeline_signature_accepts_aligned_manifest(
        self,
    ) -> None:
        self.assertIn(
            (
                "aligned_receptor_"
                "ensemble_json"
            ),
            inspect.signature(
                run_pipeline
            ).parameters,
        )


if __name__ == "__main__":
    unittest.main()
