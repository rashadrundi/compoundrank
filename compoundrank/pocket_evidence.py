from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Iterable

from .models import PocketDefinition


SCHEMA_VERSION = "pocket_evidence.v0.1"

SELECTION_MODES = {
    "report_only",
    "prioritize_supported",
}

EVIDENCE_ORIGINS = {
    "curated_annotation",
    "homolog_transfer",
    "expert_input",
    "reference_ligand_posthoc",
}


def normalize_residue_id(
    value: object,
) -> str:
    text = str(value or "").strip()

    match = re.fullmatch(
        r"\s*([A-Za-z]{3})\s*:\s*"
        r"([^:]*)\s*:\s*"
        r"(-?\d+)([A-Za-z]?)\s*",
        text,
    )

    if match is None:
        raise ValueError(
            "Residues must use RES:CHAIN:NUMBER format, "
            f"for example ARG:B:118; received {value!r}"
        )

    residue_name = match.group(1).upper()
    chain = (
        match.group(2).strip().upper()
        or "_"
    )
    residue_number = match.group(3)
    insertion_code = match.group(4).upper()

    return (
        f"{residue_name}:"
        f"{chain}:"
        f"{residue_number}"
        f"{insertion_code}"
    )


def load_pocket_evidence(
    path: Path,
) -> dict[str, Any]:
    evidence_path = Path(path)

    if not evidence_path.is_file():
        raise FileNotFoundError(
            evidence_path
        )

    payload = json.loads(
        evidence_path.read_text(
            encoding="utf-8"
        )
    )

    if not isinstance(payload, dict):
        raise ValueError(
            "Pocket evidence must be a JSON object."
        )

    nested = payload.get(
        "pocket_residue_evidence"
    )

    if isinstance(nested, dict):
        payload = nested

    schema_version = str(
        payload.get("schema_version")
        or SCHEMA_VERSION
    )

    if schema_version != SCHEMA_VERSION:
        raise ValueError(
            "Unsupported pocket-evidence schema: "
            f"{schema_version}"
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

    evidence_origin = str(
        payload.get("evidence_origin")
        or ""
    )

    if evidence_origin not in EVIDENCE_ORIGINS:
        raise ValueError(
            "evidence_origin must be one of: "
            + ", ".join(
                sorted(EVIDENCE_ORIGINS)
            )
        )

    if (
        evidence_origin
        == "reference_ligand_posthoc"
        and selection_mode
        == "prioritize_supported"
    ):
        raise ValueError(
            "Post-hoc reference-ligand evidence cannot "
            "be used for pocket selection."
        )

    raw_residues = payload.get(
        "residues"
    )

    if not isinstance(raw_residues, list):
        raise ValueError(
            "Pocket evidence requires a residues list."
        )

    residues = sorted(
        {
            normalize_residue_id(value)
            for value in raw_residues
        }
    )

    if not residues:
        raise ValueError(
            "Pocket evidence contains no usable residues."
        )

    confidence = str(
        payload.get("confidence")
        or "unspecified"
    ).lower()

    evidence_id = str(
        payload.get("evidence_id")
        or evidence_path.stem
    )

    return {
        "schema_version": SCHEMA_VERSION,
        "evidence_id": evidence_id,
        "evidence_origin": evidence_origin,
        "selection_mode": selection_mode,
        "confidence": confidence,
        "residues": residues,
        "source": payload.get("source"),
        "notes": list(
            payload.get("notes") or []
        ),
        "input_path": str(evidence_path),
    }


def _parse_pdb_residue_ids(
    path: Path,
) -> set[str]:
    residues: set[str] = set()

    for line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        if not line.startswith(
            ("ATOM", "HETATM")
        ):
            continue

        residue_name = (
            line[17:20].strip().upper()
        )
        chain = (
            line[21:22].strip()
            or "_"
        )
        residue_number = (
            line[22:26].strip()
        )
        insertion_code = (
            line[26:27].strip().upper()
        )

        if (
            not residue_name
            or not residue_number
        ):
            continue

        residues.add(
            f"{residue_name}:"
            f"{chain}:"
            f"{residue_number}"
            f"{insertion_code}"
        )

    return residues


def _fpocket_number(
    pocket_id: str,
) -> int | None:
    match = re.search(
        r"_pocket_(\d+)$",
        pocket_id,
    )

    if match is None:
        return None

    return int(match.group(1))


def collect_pocket_lining_residues(
    pockets: Iterable[PocketDefinition],
    *,
    fpocket_output_dir: Path,
) -> tuple[
    dict[str, set[str]],
    dict[str, str],
]:
    pocket_list = list(pockets)

    by_id = {
        pocket.pocket_id: pocket
        for pocket in pocket_list
    }

    cache: dict[str, set[str]] = {}
    statuses: dict[str, str] = {}

    def collect(
        pocket: PocketDefinition,
        stack: set[str],
    ) -> set[str]:
        if pocket.pocket_id in cache:
            return set(
                cache[pocket.pocket_id]
            )

        if pocket.pocket_id in stack:
            raise ValueError(
                "Circular merged-pocket definition: "
                f"{pocket.pocket_id}"
            )

        if pocket.merged_from:
            merged_residues: set[str] = set()
            next_stack = {
                *stack,
                pocket.pocket_id,
            }

            missing_components = []

            for component_id in (
                pocket.merged_from
            ):
                component = by_id.get(
                    component_id
                )

                if component is None:
                    missing_components.append(
                        component_id
                    )
                    continue

                merged_residues.update(
                    collect(
                        component,
                        next_stack,
                    )
                )

            cache[pocket.pocket_id] = (
                merged_residues
            )

            if missing_components:
                statuses[pocket.pocket_id] = (
                    "partial_missing_components:"
                    + ",".join(
                        missing_components
                    )
                )
            elif merged_residues:
                statuses[pocket.pocket_id] = (
                    "merged_component_union"
                )
            else:
                statuses[pocket.pocket_id] = (
                    "merged_components_without_residues"
                )

            return set(merged_residues)

        number = _fpocket_number(
            pocket.pocket_id
        )

        if number is None:
            cache[pocket.pocket_id] = set()
            statuses[pocket.pocket_id] = (
                "not_fpocket_derived"
            )
            return set()

        atom_path = (
            fpocket_output_dir
            / "pockets"
            / f"pocket{number}_atm.pdb"
        )

        if not atom_path.is_file():
            cache[pocket.pocket_id] = set()
            statuses[pocket.pocket_id] = (
                "lining_residue_file_missing"
            )
            return set()

        residues = _parse_pdb_residue_ids(
            atom_path
        )

        cache[pocket.pocket_id] = residues
        statuses[pocket.pocket_id] = (
            "fpocket_atom_file"
            if residues
            else "fpocket_atom_file_empty"
        )

        return set(residues)

    for pocket in pocket_list:
        collect(pocket, set())

    return cache, statuses


def score_pocket_biological_evidence(
    pockets: Iterable[PocketDefinition],
    *,
    fpocket_output_dir: Path,
    evidence: dict[str, Any],
) -> tuple[
    dict[str, dict[str, Any]],
    dict[str, Any],
]:
    pocket_list = list(pockets)

    evidence_residues = set(
        evidence["residues"]
    )

    lining_by_id, statuses = (
        collect_pocket_lining_residues(
            pocket_list,
            fpocket_output_dir=(
                fpocket_output_dir
            ),
        )
    )

    use_for_selection = (
        evidence["selection_mode"]
        == "prioritize_supported"
    )

    scores: dict[
        str,
        dict[str, Any],
    ] = {}

    report_rows = []

    for pocket in pocket_list:
        lining = lining_by_id.get(
            pocket.pocket_id,
            set(),
        )

        overlap = (
            lining & evidence_residues
        )
        union = (
            lining | evidence_residues
        )

        recall = (
            len(overlap)
            / len(evidence_residues)
            if evidence_residues
            else 0.0
        )

        precision = (
            len(overlap)
            / len(lining)
            if lining
            else 0.0
        )

        jaccard = (
            len(overlap)
            / len(union)
            if union
            else 0.0
        )

        row = {
            "biological_evidence_available": True,
            "biological_evidence_used_for_selection": (
                use_for_selection
            ),
            "biological_evidence_supported": (
                bool(overlap)
            ),
            "biological_evidence_id": (
                evidence["evidence_id"]
            ),
            "biological_evidence_origin": (
                evidence["evidence_origin"]
            ),
            "biological_evidence_confidence": (
                evidence["confidence"]
            ),
            "biological_evidence_overlap_count": (
                len(overlap)
            ),
            "biological_evidence_recall": recall,
            "biological_evidence_precision": (
                precision
            ),
            "biological_evidence_jaccard": (
                jaccard
            ),
            "biological_evidence_matched_residues": (
                sorted(overlap)
            ),
            "pocket_lining_residue_count": (
                len(lining)
            ),
            "pocket_lining_residues": (
                sorted(lining)
            ),
            "pocket_lining_source_status": (
                statuses.get(
                    pocket.pocket_id,
                    "unknown",
                )
            ),
        }

        scores[pocket.pocket_id] = row

        report_rows.append(
            {
                "pocket_id": pocket.pocket_id,
                "pocket_rank": (
                    pocket.pocket_rank
                ),
                "fpocket_score": (
                    pocket.fpocket_score
                ),
                "merged_from": list(
                    pocket.merged_from
                ),
                **row,
            }
        )

    report_rows.sort(
        key=lambda row: (
            row[
                "biological_evidence_recall"
            ],
            row[
                "biological_evidence_jaccard"
            ],
            row[
                "biological_evidence_precision"
            ],
            row[
                "biological_evidence_overlap_count"
            ],
            -int(row["pocket_rank"]),
        ),
        reverse=True,
    )

    for rank, row in enumerate(
        report_rows,
        start=1,
    ):
        row[
            "biological_evidence_rank"
        ] = rank

        scores[
            str(row["pocket_id"])
        ][
            "biological_evidence_rank"
        ] = rank

    report = {
        "schema_version": (
            "pocket_biological_evidence.v0.1"
        ),
        "evidence": evidence,
        "selection_mode": (
            evidence["selection_mode"]
        ),
        "used_for_selection": (
            use_for_selection
        ),
        "evidence_residue_count": len(
            evidence_residues
        ),
        "pocket_count": len(
            report_rows
        ),
        "pockets": report_rows,
    }

    return scores, report


def write_pocket_biological_evidence(
    path: Path,
    payload: dict[str, Any],
) -> Path:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return path
