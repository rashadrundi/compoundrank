from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .pocket_evidence import (
    SCHEMA_VERSION as POCKET_EVIDENCE_SCHEMA,
    SELECTION_MODES,
)
from .residue_transfer import (
    transfer_reference_residues_to_structure,
)


REFERENCE_SCHEMA_VERSION = (
    "functional_site_reference.v0.1"
)

TRANSFER_SCHEMA_VERSION = (
    "homolog_pocket_evidence_transfer.v0.1"
)

SELECTION_ALIGNMENT_GRADES = {
    "strong",
    "moderate",
}


def _normalize_reference_record(
    payload: dict[str, Any],
) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise ValueError(
            "Functional-site reference must be "
            "a JSON object."
        )

    schema_version = str(
        payload.get("schema_version")
        or REFERENCE_SCHEMA_VERSION
    )

    if (
        schema_version
        != REFERENCE_SCHEMA_VERSION
    ):
        raise ValueError(
            "Unsupported functional-site "
            f"reference schema: {schema_version}"
        )

    evidence_id = str(
        payload.get("evidence_id")
        or ""
    ).strip()

    if not evidence_id:
        raise ValueError(
            "Functional-site reference requires "
            "an evidence_id."
        )

    selection_mode = str(
        payload.get("selection_mode")
        or "report_only"
    )

    if selection_mode not in SELECTION_MODES:
        raise ValueError(
            "selection_mode must be one of: "
            + ", ".join(
                sorted(SELECTION_MODES)
            )
        )

    residues = payload.get("residues")

    if not isinstance(residues, list):
        raise ValueError(
            "Functional-site reference requires "
            "a residues list."
        )

    if not residues:
        raise ValueError(
            "Functional-site reference contains "
            "no residues."
        )

    normalized_residues = []

    for residue in residues:
        if not isinstance(residue, dict):
            raise ValueError(
                "Each functional-site residue "
                "must be a JSON object."
            )

        if (
            "residue_number" not in residue
            and "sequence_position"
            not in residue
            and "position" not in residue
        ):
            raise ValueError(
                "Each functional-site residue "
                "requires residue_number or "
                "sequence_position."
            )

        normalized_residues.append(
            dict(residue)
        )

    return {
        "schema_version": (
            REFERENCE_SCHEMA_VERSION
        ),
        "evidence_id": evidence_id,
        "selection_mode": selection_mode,
        "confidence": str(
            payload.get("confidence")
            or "unspecified"
        ).lower(),
        "reference_chain": (
            str(
                payload.get(
                    "reference_chain"
                )
                or ""
            ).strip().upper()
            or None
        ),
        "residues": normalized_residues,
        "source": payload.get("source"),
        "notes": list(
            payload.get("notes") or []
        ),
    }


def load_functional_site_reference(
    path: Path,
) -> dict[str, Any]:
    reference_path = Path(path)

    if not reference_path.is_file():
        raise FileNotFoundError(
            reference_path
        )

    payload = json.loads(
        reference_path.read_text(
            encoding="utf-8"
        )
    )

    normalized = (
        _normalize_reference_record(
            payload
        )
    )

    normalized["input_path"] = str(
        reference_path
    )

    return normalized


def _effective_confidence(
    *,
    requested_confidence: str,
    homolog_grade: str,
    mapping_fraction: float,
    mapped_count: int,
) -> str:
    if (
        homolog_grade == "strong"
        and mapping_fraction >= 0.80
        and mapped_count >= 3
    ):
        return (
            requested_confidence
            if requested_confidence
            in {"high", "moderate"}
            else "high"
        )

    if (
        homolog_grade
        in {"strong", "moderate"}
        and mapping_fraction >= 0.50
        and mapped_count >= 2
    ):
        return "moderate"

    return "low"


def build_homolog_pocket_evidence(
    *,
    reference_record: dict[str, Any],
    reference_sequence: str,
    submitted_sequence: str,
    reference_pdb: Path,
    receptor_pdb: Path,
    reference_chain_id: str | None = None,
    receptor_chain_id: str | None = None,
    minimum_mapping_fraction: float = 0.50,
    minimum_mapped_residues: int = 2,
    minimum_structure_identity: float = 0.90,
    minimum_structure_coverage: float = 0.80,
) -> dict[str, Any]:
    """Transfer homolog functional residues into pocket evidence.

    Selection is enabled only when the requested mode is
    ``prioritize_supported`` and all conservative transfer
    thresholds are satisfied. Otherwise evidence is retained
    in ``report_only`` mode.
    """

    record = _normalize_reference_record(
        reference_record
    )

    selected_reference_chain = (
        reference_chain_id
        or record["reference_chain"]
    )

    transfer = (
        transfer_reference_residues_to_structure(
            reference_sequence=(
                reference_sequence
            ),
            submitted_sequence=(
                submitted_sequence
            ),
            reference_pdb=reference_pdb,
            reference_chain_id=(
                selected_reference_chain
            ),
            receptor_pdb=receptor_pdb,
            chain_id=receptor_chain_id,
            reference_residues=(
                record["residues"]
            ),
        )
    )

    mapped_count = int(
        transfer[
            "mapped_conserved_count"
        ]
    )

    mapping_fraction = float(
        transfer["mapping_fraction"]
    )

    homolog_alignment = transfer[
        "reference_to_submitted_alignment"
    ]

    structure_alignment = transfer[
        "submitted_to_structure_alignment"
    ]

    homolog_grade = str(
        homolog_alignment.get(
            "alignment_grade"
        )
        or "exploratory"
    )

    structure_identity = float(
        structure_alignment.get(
            "identity"
        )
        or 0.0
    )

    structure_coverage = float(
        structure_alignment.get(
            "reference_coverage"
        )
        or 0.0
    )

    checks = {
        "selection_requested": (
            record["selection_mode"]
            == "prioritize_supported"
        ),
        "homolog_alignment_supported": (
            homolog_grade
            in SELECTION_ALIGNMENT_GRADES
        ),
        "mapping_fraction_supported": (
            mapping_fraction
            >= minimum_mapping_fraction
        ),
        "mapped_residue_count_supported": (
            mapped_count
            >= minimum_mapped_residues
        ),
        "structure_identity_supported": (
            structure_identity
            >= minimum_structure_identity
        ),
        "structure_coverage_supported": (
            structure_coverage
            >= minimum_structure_coverage
        ),
    }

    use_for_selection = all(
        checks.values()
    )

    effective_selection_mode = (
        "prioritize_supported"
        if use_for_selection
        else "report_only"
    )

    effective_confidence = (
        _effective_confidence(
            requested_confidence=(
                record["confidence"]
            ),
            homolog_grade=(
                homolog_grade
            ),
            mapping_fraction=(
                mapping_fraction
            ),
            mapped_count=mapped_count,
        )
    )

    notes = list(record["notes"])

    notes.append(
        "Functional residues were transferred "
        "through reference numbering, homolog "
        "sequence alignment, and receptor "
        "structure alignment."
    )

    if (
        record["selection_mode"]
        == "prioritize_supported"
        and not use_for_selection
    ):
        failed_checks = sorted(
            name
            for name, passed in checks.items()
            if not passed
            and name
            != "selection_requested"
        )

        notes.append(
            "Requested biological prioritization "
            "was downgraded to report_only "
            "because these checks failed: "
            + ", ".join(failed_checks)
        )

    return {
        "schema_version": (
            POCKET_EVIDENCE_SCHEMA
        ),
        "evidence_id": (
            record["evidence_id"]
        ),
        "evidence_origin": (
            "homolog_transfer"
        ),
        "selection_mode": (
            effective_selection_mode
        ),
        "confidence": (
            effective_confidence
        ),
        "residues": list(
            transfer["residues"]
        ),
        "source": {
            "functional_site_source": (
                record["source"]
            ),
            "reference_pdb": str(
                Path(reference_pdb)
            ),
            "reference_chain": (
                transfer[
                    "reference_numbering_map"
                ]["selected_chain"]
                if transfer.get(
                    "reference_numbering_map"
                )
                else None
            ),
            "receptor_pdb": str(
                Path(receptor_pdb)
            ),
            "receptor_chain": (
                transfer[
                    "selected_structure_chain"
                ]
            ),
            "homolog_alignment_grade": (
                homolog_grade
            ),
            "homolog_identity": (
                homolog_alignment[
                    "identity"
                ]
            ),
            "homolog_reference_coverage": (
                homolog_alignment[
                    "reference_coverage"
                ]
            ),
            "structure_identity": (
                structure_identity
            ),
            "structure_submitted_coverage": (
                structure_coverage
            ),
        },
        "notes": notes,
        "transfer_summary": {
            "schema_version": (
                TRANSFER_SCHEMA_VERSION
            ),
            "requested_selection_mode": (
                record["selection_mode"]
            ),
            "effective_selection_mode": (
                effective_selection_mode
            ),
            "selection_checks": checks,
            "thresholds": {
                "minimum_mapping_fraction": (
                    minimum_mapping_fraction
                ),
                "minimum_mapped_residues": (
                    minimum_mapped_residues
                ),
                "minimum_structure_identity": (
                    minimum_structure_identity
                ),
                "minimum_structure_coverage": (
                    minimum_structure_coverage
                ),
            },
            "requested_residue_count": (
                transfer[
                    "requested_residue_count"
                ]
            ),
            "mapped_conserved_count": (
                mapped_count
            ),
            "mapping_fraction": (
                mapping_fraction
            ),
            "status_counts": (
                transfer["status_counts"]
            ),
        },
        "transfer_report": transfer,
    }


def write_homolog_pocket_evidence(
    output_path: Path,
    *,
    reference_record: dict[str, Any],
    reference_sequence: str,
    submitted_sequence: str,
    reference_pdb: Path,
    receptor_pdb: Path,
    reference_chain_id: str | None = None,
    receptor_chain_id: str | None = None,
) -> dict[str, Any]:
    output = (
        build_homolog_pocket_evidence(
            reference_record=(
                reference_record
            ),
            reference_sequence=(
                reference_sequence
            ),
            submitted_sequence=(
                submitted_sequence
            ),
            reference_pdb=reference_pdb,
            receptor_pdb=receptor_pdb,
            reference_chain_id=(
                reference_chain_id
            ),
            receptor_chain_id=(
                receptor_chain_id
            ),
        )
    )

    path = Path(output_path)

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        json.dumps(
            output,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return output
