from __future__ import annotations

import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any


SOURCE_SCHEMA_VERSION = (
    "aligned_receptor_ensemble.v0.1"
)

AUDIT_SCHEMA_VERSION = (
    "aligned_receptor_ensemble_input.v0.1"
)

ALIGNMENT_METHOD = (
    "chain_residue_ca_kabsch.v1"
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
            f"{label} must be an object"
        )

    return dict(value)


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
            "nonempty string"
        )

    return value.strip()


def _required_integer(
    mapping: dict[str, Any],
    key: str,
    label: str,
) -> int:
    value = mapping.get(key)

    if (
        not isinstance(value, int)
        or isinstance(value, bool)
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
        not isinstance(value, (int, float))
        or isinstance(value, bool)
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
    value: str,
    *,
    label: str,
    expected_checksum: str,
    verify_checksum: bool,
) -> Path:
    path = (
        Path(value)
        .expanduser()
        .resolve()
    )

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):
        raise FileNotFoundError(
            f"{label} is missing or empty: "
            f"{path}"
        )

    if verify_checksum:
        actual_checksum = _sha256(path)

        if actual_checksum != expected_checksum:
            raise ValueError(
                f"{label} checksum mismatch: "
                f"expected {expected_checksum}; "
                f"found {actual_checksum}"
            )

    return path


def _validate_conformer_id(
    conformer_id: str,
) -> None:
    if conformer_id == "submitted_receptor":
        raise ValueError(
            "Snapshot conformer ID cannot use "
            "the reserved submitted_receptor ID"
        )

    if any(
        separator in conformer_id
        for separator in (
            "/",
            "\\",
        )
    ):
        raise ValueError(
            "Snapshot conformer ID cannot "
            "contain path separators"
        )


def validate_receptor_ensemble_options(
    report_only_manifest: Path | None,
    aligned_manifest: Path | None,
) -> None:
    if (
        report_only_manifest is not None
        and aligned_manifest is not None
    ):
        raise ValueError(
            "Use either the report-only receptor "
            "ensemble manifest or the aligned "
            "receptor ensemble manifest, not both"
        )


def load_aligned_receptor_ensemble_manifest(
    manifest_path: Path,
    *,
    submitted_receptor_pdb: Path,
    verify_checksums: bool = True,
) -> dict[str, Any]:
    path = (
        Path(manifest_path)
        .expanduser()
        .resolve()
    )

    submitted_receptor = (
        Path(submitted_receptor_pdb)
        .expanduser()
        .resolve()
    )

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):
        raise FileNotFoundError(
            "Aligned receptor ensemble manifest "
            f"is missing or empty: {path}"
        )

    if (
        not submitted_receptor.is_file()
        or submitted_receptor.stat().st_size == 0
    ):
        raise FileNotFoundError(
            "Submitted receptor is missing or "
            f"empty: {submitted_receptor}"
        )

    try:
        payload = json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )
    except json.JSONDecodeError as error:
        raise ValueError(
            "Aligned receptor ensemble manifest "
            f"is not valid JSON: {path}"
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

    if schema_version != SOURCE_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported aligned receptor "
            f"ensemble schema: {schema_version!r}; "
            f"expected {SOURCE_SCHEMA_VERSION!r}"
        )

    status = _required_string(
        manifest,
        "status",
        "manifest",
    )

    if status != "complete":
        raise ValueError(
            "Aligned receptor ensemble must "
            "have status 'complete'"
        )

    selection_mode = _required_string(
        manifest,
        "selection_mode",
        "manifest",
    )

    if selection_mode != "report_only":
        raise ValueError(
            "Aligned receptor ensemble must "
            "retain report_only selection mode"
        )

    alignment_method = _required_string(
        manifest,
        "alignment_method",
        "manifest",
    )

    if alignment_method != ALIGNMENT_METHOD:
        raise ValueError(
            "Unsupported aligned receptor "
            f"alignment method: {alignment_method}"
        )

    docking_behavior = _required_string(
        manifest,
        "docking_behavior",
        "manifest",
    )

    if docking_behavior != "not_enabled":
        raise ValueError(
            "Aligned generation artifact must "
            "have docking_behavior 'not_enabled'"
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
            "aligned_reference_path",
            "manifest.reference",
        ),
        label=(
            "manifest.reference."
            "aligned_reference_path"
        ),
        expected_checksum=(
            reference_checksum
        ),
        verify_checksum=(
            verify_checksums
        ),
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

    submitted_checksum = _sha256(
        submitted_receptor
    )

    if submitted_checksum != reference_checksum:
        raise ValueError(
            "The submitted receptor does not "
            "exactly match the aligned ensemble "
            "reference. Shared pocket coordinates "
            "cannot be reused safely."
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

    if snapshot_count <= 0:
        raise ValueError(
            "Aligned receptor ensemble contains "
            "no snapshots"
        )

    normalized_snapshots: list[
        dict[str, Any]
    ] = []

    seen_ids: set[str] = set()

    for index, raw_snapshot in enumerate(
        snapshots_value
    ):
        label = (
            f"manifest.snapshots[{index}]"
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

        _validate_conformer_id(
            snapshot_id
        )

        if snapshot_id in seen_ids:
            raise ValueError(
                "Duplicate aligned receptor "
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

        if snapshot_status != "aligned":
            raise ValueError(
                f"{label}.status must be 'aligned'"
            )

        checksum = _required_string(
            snapshot,
            "aligned_checksum_sha256",
            label,
        )

        aligned_path = _validated_file(
            _required_string(
                snapshot,
                "aligned_path",
                label,
            ),
            label=f"{label}.aligned_path",
            expected_checksum=checksum,
            verify_checksum=(
                verify_checksums
            ),
        )

        matched_ca_atoms = _required_integer(
            snapshot,
            "matched_ca_atoms",
            label,
        )

        if matched_ca_atoms <= 0:
            raise ValueError(
                f"{label}.matched_ca_atoms must "
                "be greater than zero"
            )

        coverage = _required_finite_float(
            snapshot,
            "reference_coverage_fraction",
            label,
        )

        if not 0.0 < coverage <= 1.0:
            raise ValueError(
                f"{label}.reference_coverage_"
                "fraction must be greater than "
                "zero and no greater than one"
            )

        kabsch_rmsd = (
            _required_finite_float(
                snapshot,
                "kabsch_ca_rmsd_angstrom",
                label,
            )
        )

        raw_aligned_rmsd = (
            _required_finite_float(
                snapshot,
                (
                    "raw_ca_rmsd_after_"
                    "alignment_angstrom"
                ),
                label,
            )
        )

        if (
            kabsch_rmsd < 0.0
            or raw_aligned_rmsd < 0.0
        ):
            raise ValueError(
                f"{label} RMSD values cannot "
                "be negative"
            )

        if abs(
            raw_aligned_rmsd
            - kabsch_rmsd
        ) > 0.005:
            raise ValueError(
                f"{label} does not appear to be "
                "stored in the reference "
                "coordinate frame"
            )

        normalized_snapshots.append(
            {
                "conformer_id": snapshot_id,
                "status": "accepted",
                "aligned_path": str(
                    aligned_path
                ),
                "checksum_sha256": checksum,
                "matched_ca_atoms": (
                    matched_ca_atoms
                ),
                (
                    "reference_coverage_"
                    "fraction"
                ): coverage,
                (
                    "kabsch_ca_rmsd_"
                    "angstrom"
                ): kabsch_rmsd,
                (
                    "raw_ca_rmsd_after_"
                    "alignment_angstrom"
                ): raw_aligned_rmsd,
            }
        )

    return {
        "schema_version": (
            AUDIT_SCHEMA_VERSION
        ),
        "status": "accepted",
        "selection_mode": "report_only",
        "docking_behavior": (
            "validated_not_prepared"
        ),
        "source_manifest": str(
            path
        ),
        "source_manifest_checksum_sha256": (
            _sha256(path)
        ),
        "alignment_method": (
            alignment_method
        ),
        "reference": {
            "submitted_receptor_path": str(
                submitted_receptor
            ),
            "aligned_reference_path": str(
                reference_path
            ),
            "checksum_sha256": (
                reference_checksum
            ),
            "ca_atoms": (
                reference_ca_atoms
            ),
        },
        "snapshot_count": len(
            normalized_snapshots
        ),
        "snapshots": (
            normalized_snapshots
        ),
        "limitations": [
            (
                "This stage validates aligned "
                "receptor conformers only."
            ),
            (
                "The receptor conformers have not "
                "yet been prepared for GNINA."
            ),
            (
                "Docking behavior remains limited "
                "to the submitted receptor."
            ),
        ],
    }


def record_aligned_receptor_ensemble_input(
    *,
    manifest_path: Path,
    submitted_receptor_pdb: Path,
    output_dir: Path,
    verify_checksums: bool = True,
    overwrite: bool = False,
) -> tuple[dict[str, Any], Path]:
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
            "Aligned receptor ensemble audit "
            f"directory is not empty: {destination}"
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

    audit = (
        load_aligned_receptor_ensemble_manifest(
            manifest_path,
            submitted_receptor_pdb=(
                submitted_receptor_pdb
            ),
            verify_checksums=(
                verify_checksums
            ),
        )
    )

    source_copy = (
        destination
        / "aligned_receptor_ensemble.json"
    )

    shutil.copy2(
        Path(
            audit["source_manifest"]
        ),
        source_copy,
    )

    audit_path = (
        destination
        / "aligned_receptor_ensemble_input.json"
    )

    audit["outputs"] = {
        "source_manifest_copy": str(
            source_copy
        ),
        "audit_json": str(
            audit_path
        ),
    }

    audit_path.write_text(
        json.dumps(
            audit,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return audit, audit_path
