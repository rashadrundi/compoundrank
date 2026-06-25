from __future__ import annotations

import hashlib
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

from .models import PreparedReceptor
from .receptor import prepare_receptor


SCHEMA_VERSION = (
    "receptor_conformer_preparation.v0.1"
)

PrepareReceptorFunction = Callable[
    ...,
    PreparedReceptor,
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def prepare_receptor_conformers(
    *,
    submitted_receptor_pdb: Path,
    aligned_ensemble: dict[str, Any],
    cache_root: Path,
    ph: float,
    pdb2pqr_bin: str,
    meeko_receptor_bin: str,
    prepare_receptor_fn: (
        PrepareReceptorFunction
    ) = prepare_receptor,
) -> list[
    tuple[str, PreparedReceptor]
]:
    submitted_path = (
        Path(submitted_receptor_pdb)
        .expanduser()
        .resolve()
    )

    conformer_sources: list[
        tuple[str, Path]
    ] = [
        (
            "submitted_receptor",
            submitted_path,
        )
    ]

    snapshots = aligned_ensemble.get(
        "snapshots"
    )

    if not isinstance(
        snapshots,
        list,
    ):
        raise ValueError(
            "Aligned receptor ensemble audit "
            "must contain a snapshots list"
        )

    seen_ids = {
        "submitted_receptor"
    }

    for index, snapshot_value in enumerate(
        snapshots
    ):
        if not isinstance(
            snapshot_value,
            dict,
        ):
            raise ValueError(
                "Aligned receptor ensemble "
                f"snapshot {index} must be an object"
            )

        conformer_id_value = (
            snapshot_value.get(
                "conformer_id"
            )
        )

        if (
            not isinstance(
                conformer_id_value,
                str,
            )
            or not conformer_id_value.strip()
        ):
            raise ValueError(
                "Aligned receptor snapshot "
                f"{index} has no conformer ID"
            )

        conformer_id = (
            conformer_id_value.strip()
        )

        if conformer_id in seen_ids:
            raise ValueError(
                "Duplicate receptor conformer ID: "
                f"{conformer_id}"
            )

        seen_ids.add(
            conformer_id
        )

        aligned_path_value = (
            snapshot_value.get(
                "aligned_path"
            )
        )

        if (
            not isinstance(
                aligned_path_value,
                str,
            )
            or not aligned_path_value.strip()
        ):
            raise ValueError(
                "Aligned receptor snapshot "
                f"{conformer_id} has no path"
            )

        aligned_path = (
            Path(
                aligned_path_value
            )
            .expanduser()
            .resolve()
        )

        if (
            not aligned_path.is_file()
            or aligned_path.stat().st_size
            == 0
        ):
            raise FileNotFoundError(
                "Aligned receptor snapshot is "
                f"missing or empty: {aligned_path}"
            )

        conformer_sources.append(
            (
                conformer_id,
                aligned_path,
            )
        )

    prepared_conformers: list[
        tuple[str, PreparedReceptor]
    ] = []

    for conformer_id, source_path in (
        conformer_sources
    ):
        prepared = prepare_receptor_fn(
            source_path,
            cache_root,
            ph=ph,
            pdb2pqr_bin=pdb2pqr_bin,
            meeko_receptor_bin=(
                meeko_receptor_bin
            ),
        )

        prepared_conformers.append(
            (
                conformer_id,
                prepared,
            )
        )

    return prepared_conformers


def write_receptor_conformer_preparation(
    output_path: Path,
    prepared_conformers: list[
        tuple[str, PreparedReceptor]
    ],
    *,
    aligned_ensemble: dict[str, Any],
) -> Path:
    destination = (
        Path(output_path)
        .expanduser()
        .resolve()
    )

    if not prepared_conformers:
        raise ValueError(
            "No prepared receptor conformers "
            "were supplied"
        )

    rows: list[
        dict[str, Any]
    ] = []

    seen_ids: set[str] = set()

    for conformer_id, receptor in (
        prepared_conformers
    ):
        if conformer_id in seen_ids:
            raise ValueError(
                "Duplicate prepared receptor "
                f"conformer ID: {conformer_id}"
            )

        seen_ids.add(
            conformer_id
        )

        paths = {
            "source_pdb": Path(
                receptor.source_pdb
            ),
            "prepared_pdbqt": Path(
                receptor.prepared_pdbqt
            ),
            "display_pdb": Path(
                receptor.display_pdb
            ),
        }

        for label, path in paths.items():
            if (
                not path.is_file()
                or path.stat().st_size == 0
            ):
                raise FileNotFoundError(
                    "Prepared receptor conformer "
                    f"{label} is missing or empty: "
                    f"{path}"
                )

        rows.append(
            {
                "conformer_id": (
                    conformer_id
                ),
                "role": (
                    "submitted_reference"
                    if conformer_id
                    == "submitted_receptor"
                    else "aligned_snapshot"
                ),
                "source_pdb": str(
                    paths["source_pdb"]
                    .resolve()
                ),
                "source_checksum_sha256": (
                    _sha256(
                        paths["source_pdb"]
                    )
                ),
                "prepared_pdbqt": str(
                    paths[
                        "prepared_pdbqt"
                    ].resolve()
                ),
                "display_pdb": str(
                    paths[
                        "display_pdb"
                    ].resolve()
                ),
                "cache_key": (
                    receptor.cache_key
                ),
            }
        )

    expected_count = (
        1
        + int(
            aligned_ensemble[
                "snapshot_count"
            ]
        )
    )

    if len(rows) != expected_count:
        raise ValueError(
            "Prepared receptor conformer count "
            "does not match aligned ensemble: "
            f"expected {expected_count}; "
            f"found {len(rows)}"
        )

    payload = {
        "schema_version": (
            SCHEMA_VERSION
        ),
        "status": "complete",
        "selection_mode": (
            "shared_aligned_pocket_frame"
        ),
        "docking_behavior": (
            "submitted_receptor_only"
        ),
        "source_aligned_ensemble": (
            aligned_ensemble.get(
                "source_manifest"
            )
        ),
        "prepared_conformer_count": (
            len(rows)
        ),
        "conformers": rows,
        "limitations": [
            (
                "Every listed receptor conformer "
                "has been prepared for GNINA."
            ),
            (
                "GNINA still docks only against "
                "the submitted receptor in this "
                "preparation milestone."
            ),
            (
                "No thermodynamic or population "
                "weight is assigned to any "
                "receptor conformer."
            ),
        ],
    }

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    destination.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return destination
