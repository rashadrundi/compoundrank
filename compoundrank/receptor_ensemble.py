from __future__ import annotations

import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "structure_ensemble.v0.1"
AUDIT_SCHEMA_VERSION = (
    "receptor_ensemble_input.v0.1"
)


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


def _required_mapping(
    value: object,
    label: str,
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(
            f"{label} must be a JSON object"
        )

    return value


def _required_string(
    mapping: dict[str, Any],
    key: str,
    label: str,
) -> str:
    value = mapping.get(key)

    if (
        not isinstance(value, str)
        or not value.strip()
    ):
        raise ValueError(
            f"{label}.{key} must be a "
            "non-empty string"
        )

    return value.strip()


def _required_integer(
    mapping: dict[str, Any],
    key: str,
    label: str,
) -> int:
    value = mapping.get(key)

    if (
        isinstance(value, bool)
        or not isinstance(value, int)
    ):
        raise ValueError(
            f"{label}.{key} must be an integer"
        )

    return value


def _required_finite_float(
    mapping: dict[str, Any],
    key: str,
    label: str,
) -> float:
    value = mapping.get(key)

    if (
        isinstance(value, bool)
        or not isinstance(
            value,
            (
                int,
                float,
            ),
        )
    ):
        raise ValueError(
            f"{label}.{key} must be numeric"
        )

    result = float(value)

    if not math.isfinite(result):
        raise ValueError(
            f"{label}.{key} must be finite"
        )

    return result


def _validated_file(
    raw_path: str,
    *,
    label: str,
    expected_checksum: str,
    verify_checksum: bool,
) -> Path:
    path = (
        Path(raw_path)
        .expanduser()
        .resolve()
    )

    if not path.is_absolute():
        raise ValueError(
            f"{label} must be an absolute path"
        )

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):
        raise FileNotFoundError(
            f"{label} is missing or empty: {path}"
        )

    if path.suffix.lower() != ".pdb":
        raise ValueError(
            f"{label} must reference a PDB file: "
            f"{path}"
        )

    if verify_checksum:
        actual = _sha256(path)

        if actual != expected_checksum:
            raise ValueError(
                f"{label} checksum mismatch: "
                f"expected {expected_checksum}, "
                f"found {actual}"
            )

    return path


def load_receptor_ensemble_manifest(
    manifest_path: Path,
    *,
    verify_checksums: bool = True,
) -> dict[str, Any]:
    path = (
        Path(manifest_path)
        .expanduser()
        .resolve()
    )

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):
        raise FileNotFoundError(
            f"Receptor ensemble manifest is "
            f"missing or empty: {path}"
        )

    try:
        payload = json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )
    except json.JSONDecodeError as error:
        raise ValueError(
            "Receptor ensemble manifest is "
            f"not valid JSON: {path}"
        ) from error

    manifest = _required_mapping(
        payload,
        "manifest",
    )

    schema_version = _required_string(
        manifest,
        "schema_version",
        "manifest",
    )

    if schema_version != SCHEMA_VERSION:
        raise ValueError(
            "Unsupported receptor ensemble "
            f"schema: {schema_version!r}; "
            f"expected {SCHEMA_VERSION!r}"
        )

    status = _required_string(
        manifest,
        "status",
        "manifest",
    )

    if status != "complete":
        raise ValueError(
            "Receptor ensemble manifest must "
            f"have status 'complete', found "
            f"{status!r}"
        )

    selection_mode = _required_string(
        manifest,
        "selection_mode",
        "manifest",
    )

    if selection_mode != "report_only":
        raise ValueError(
            "Receptor ensemble import currently "
            "supports report_only manifests only"
        )

    source_engine = _required_string(
        manifest,
        "source_engine",
        "manifest",
    )

    reference = _required_mapping(
        manifest.get("reference"),
        "manifest.reference",
    )

    reference_checksum = _required_string(
        reference,
        "checksum_sha256",
        "manifest.reference",
    )

    reference_path = _validated_file(
        _required_string(
            reference,
            "stored_path",
            "manifest.reference",
        ),
        label=(
            "manifest.reference.stored_path"
        ),
        expected_checksum=(
            reference_checksum
        ),
        verify_checksum=verify_checksums,
    )

    reference_ca_atoms = _required_integer(
        reference,
        "ca_atoms",
        "manifest.reference",
    )

    if reference_ca_atoms <= 0:
        raise ValueError(
            "manifest.reference.ca_atoms must "
            "be greater than zero"
        )

    snapshots_value = manifest.get(
        "snapshots"
    )

    if not isinstance(
        snapshots_value,
        list,
    ):
        raise ValueError(
            "manifest.snapshots must be a list"
        )

    snapshot_count = _required_integer(
        manifest,
        "snapshot_count",
        "manifest",
    )

    if snapshot_count != len(
        snapshots_value
    ):
        raise ValueError(
            "manifest.snapshot_count does not "
            "match the snapshots list"
        )

    accepted_expected = _required_integer(
        manifest,
        "accepted_snapshot_count",
        "manifest",
    )

    rejected_expected = _required_integer(
        manifest,
        "rejected_snapshot_count",
        "manifest",
    )

    accepted: list[
        dict[str, Any]
    ] = []

    rejected: list[
        dict[str, Any]
    ] = []

    seen_ids: set[str] = set()

    for index, raw_snapshot in enumerate(
        snapshots_value,
        start=1,
    ):
        label = (
            f"manifest.snapshots[{index - 1}]"
        )

        snapshot = _required_mapping(
            raw_snapshot,
            label,
        )

        snapshot_id = _required_string(
            snapshot,
            "snapshot_id",
            label,
        )

        if snapshot_id in seen_ids:
            raise ValueError(
                "Duplicate receptor ensemble "
                f"snapshot ID: {snapshot_id}"
            )

        seen_ids.add(
            snapshot_id
        )

        snapshot_status = _required_string(
            snapshot,
            "status",
            label,
        )

        if snapshot_status == "rejected":
            rejected.append(
                {
                    "snapshot_id": snapshot_id,
                    "status": "rejected",
                    "rejection_reason": str(
                        snapshot.get(
                            "rejection_reason",
                            "",
                        )
                    ),
                }
            )
            continue

        if snapshot_status != "accepted":
            raise ValueError(
                f"{label}.status must be "
                "'accepted' or 'rejected'"
            )

        checksum = _required_string(
            snapshot,
            "checksum_sha256",
            label,
        )

        stored_path = _validated_file(
            _required_string(
                snapshot,
                "stored_path",
                label,
            ),
            label=f"{label}.stored_path",
            expected_checksum=checksum,
            verify_checksum=(
                verify_checksums
            ),
        )

        matched_ca_atoms = (
            _required_integer(
                snapshot,
                "matched_ca_atoms",
                label,
            )
        )

        snapshot_ca_atoms = (
            _required_integer(
                snapshot,
                "snapshot_ca_atoms",
                label,
            )
        )

        snapshot_reference_atoms = (
            _required_integer(
                snapshot,
                "reference_ca_atoms",
                label,
            )
        )

        coverage = (
            _required_finite_float(
                snapshot,
                (
                    "reference_coverage_"
                    "fraction"
                ),
                label,
            )
        )

        ca_rmsd = (
            _required_finite_float(
                snapshot,
                "ca_rmsd_angstrom",
                label,
            )
        )

        if matched_ca_atoms <= 0:
            raise ValueError(
                f"{label}.matched_ca_atoms "
                "must be greater than zero"
            )

        if snapshot_ca_atoms <= 0:
            raise ValueError(
                f"{label}.snapshot_ca_atoms "
                "must be greater than zero"
            )

        if (
            snapshot_reference_atoms
            != reference_ca_atoms
        ):
            raise ValueError(
                f"{label}.reference_ca_atoms "
                "does not match the manifest "
                "reference"
            )

        if not 0.0 <= coverage <= 1.0:
            raise ValueError(
                f"{label}.reference_coverage_"
                "fraction must be between "
                "zero and one"
            )

        if ca_rmsd < 0.0:
            raise ValueError(
                f"{label}.ca_rmsd_angstrom "
                "cannot be negative"
            )

        accepted.append(
            {
                "snapshot_id": snapshot_id,
                "status": "accepted",
                "stored_path": str(
                    stored_path
                ),
                "source_path": str(
                    snapshot.get(
                        "source_path",
                        "",
                    )
                ),
                "checksum_sha256": checksum,
                "matched_ca_atoms": (
                    matched_ca_atoms
                ),
                "snapshot_ca_atoms": (
                    snapshot_ca_atoms
                ),
                "reference_ca_atoms": (
                    snapshot_reference_atoms
                ),
                (
                    "reference_coverage_"
                    "fraction"
                ): coverage,
                "ca_rmsd_angstrom": ca_rmsd,
            }
        )

    if len(accepted) != accepted_expected:
        raise ValueError(
            "manifest.accepted_snapshot_count "
            "does not match accepted records"
        )

    if len(rejected) != rejected_expected:
        raise ValueError(
            "manifest.rejected_snapshot_count "
            "does not match rejected records"
        )

    if not accepted:
        raise ValueError(
            "Receptor ensemble contains no "
            "accepted snapshots"
        )

    return {
        "schema_version": (
            AUDIT_SCHEMA_VERSION
        ),
        "status": "accepted",
        "selection_mode": "report_only",
        "docking_behavior": (
            "submitted_receptor_only"
        ),
        "source_manifest": str(
            path
        ),
        "source_manifest_checksum_sha256": (
            _sha256(path)
        ),
        "source_engine": source_engine,
        "reference": {
            "stored_path": str(
                reference_path
            ),
            "checksum_sha256": (
                reference_checksum
            ),
            "ca_atoms": (
                reference_ca_atoms
            ),
        },
        "snapshot_count": (
            snapshot_count
        ),
        "accepted_snapshot_count": (
            len(accepted)
        ),
        "rejected_snapshot_count": (
            len(rejected)
        ),
        "accepted_snapshots": accepted,
        "rejected_snapshots": rejected,
        "limitations": [
            (
                "The ensemble is imported for "
                "validation and audit only."
            ),
            (
                "Docking continues to use the "
                "submitted receptor until "
                "multi-conformer docking is "
                "implemented separately."
            ),
        ],
    }


def record_receptor_ensemble_input(
    manifest_path: Path,
    output_dir: Path,
    *,
    verify_checksums: bool = True,
    overwrite: bool = False,
) -> Path:
    source = (
        Path(manifest_path)
        .expanduser()
        .resolve()
    )

    destination = (
        Path(output_dir)
        .expanduser()
        .resolve()
    )

    if (
        destination.exists()
        and any(destination.iterdir())
        and not overwrite
    ):
        raise FileExistsError(
            "Receptor ensemble audit directory "
            f"is not empty: {destination}"
        )

    if (
        overwrite
        and destination.exists()
    ):
        shutil.rmtree(
            destination
        )

    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    audit = load_receptor_ensemble_manifest(
        source,
        verify_checksums=(
            verify_checksums
        ),
    )

    manifest_copy = (
        destination
        / "source_structure_ensemble.json"
    )

    shutil.copy2(
        source,
        manifest_copy,
    )

    audit[
        "copied_source_manifest"
    ] = str(
        manifest_copy
    )

    output_path = (
        destination
        / "receptor_ensemble_input.json"
    )

    output_path.write_text(
        json.dumps(
            audit,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return output_path
