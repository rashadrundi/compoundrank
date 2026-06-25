from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable

from .models import PocketDefinition, PoseRecord


def _number_or_none(
    value: object,
) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)

    return None


def _pose_selection_key(
    record: PoseRecord,
) -> tuple[float, float, float]:
    cnn_score = _number_or_none(
        record.cnn_score
    )
    cnn_affinity = _number_or_none(
        record.cnn_affinity
    )
    minimized_affinity = _number_or_none(
        record.minimized_affinity
    )

    return (
        (
            cnn_score
            if cnn_score is not None
            else float("-inf")
        ),
        (
            cnn_affinity
            if cnn_affinity is not None
            else float("-inf")
        ),
        (
            -minimized_affinity
            if minimized_affinity is not None
            else float("-inf")
        ),
    )


def best_pose_record(
    records: Iterable[PoseRecord],
) -> PoseRecord | None:
    record_list = list(records)

    if not record_list:
        return None

    return max(
        record_list,
        key=_pose_selection_key,
    )


def summarize_pocket_attempt(
    *,
    ligand_name: str,
    receptor_conformer_id: str = (
        "submitted_receptor"
    ),
    pocket: PocketDefinition,
    raw_records: Iterable[PoseRecord],
    accepted_records: Iterable[PoseRecord],
    rejected_pose_count: int,
) -> dict[str, Any]:
    raw_list = list(raw_records)
    accepted_list = list(
        accepted_records
    )

    if accepted_list:
        scoring_records = accepted_list
        score_source = "accepted_poses"
    elif raw_list:
        scoring_records = raw_list
        score_source = "raw_pose_fallback"
    else:
        scoring_records = []
        score_source = "none"

    best_record = best_pose_record(
        scoring_records
    )

    return {
        "compound": ligand_name,
        "receptor_conformer_id": (
            receptor_conformer_id
        ),
        "pocket_id": pocket.pocket_id,
        "pocket_rank": pocket.pocket_rank,
        "fpocket_score": (
            pocket.fpocket_score
        ),
        "pocket_source": (
            pocket.source or pocket.mode
        ),
        "raw_poses": len(raw_list),
        "accepted_poses": len(
            accepted_list
        ),
        "rejected_poses": int(
            rejected_pose_count
        ),
        "score_source": score_source,
        "top_cnn_score": (
            best_record.cnn_score
            if best_record is not None
            else None
        ),
        "top_cnn_affinity": (
            best_record.cnn_affinity
            if best_record is not None
            else None
        ),
        "top_minimized_affinity": (
            best_record.minimized_affinity
            if best_record is not None
            else None
        ),
        "top_seed": (
            best_record.seed
            if best_record is not None
            else None
        ),
        "top_pose_number": (
            best_record.pose_number
            if best_record is not None
            else None
        ),
        "selection_rank": None,
        "selected": False,
    }


def _pose_validity_priority(
    row: dict[str, Any],
) -> float:
    score_source = str(
        row.get("score_source") or ""
    )

    accepted_count = _number_or_none(
        row.get("accepted_poses")
    )

    raw_count = _number_or_none(
        row.get("raw_poses")
    )

    if (
        score_source == "accepted_poses"
        or (
            accepted_count is not None
            and accepted_count > 0
        )
    ):
        return 2.0

    if (
        score_source == "raw_pose_fallback"
        or (
            raw_count is not None
            and raw_count > 0
        )
    ):
        return 1.0

    return 0.0


def _row_selection_key(
    row: dict[str, Any],
) -> tuple[
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
]:
    cnn_score = _number_or_none(
        row.get("top_cnn_score")
    )
    cnn_affinity = _number_or_none(
        row.get("top_cnn_affinity")
    )
    minimized_affinity = _number_or_none(
        row.get(
            "top_minimized_affinity"
        )
    )
    pocket_rank = _number_or_none(
        row.get("pocket_rank")
    )

    evidence_used = bool(
        row.get(
            "biological_evidence_used_for_selection"
        )
    )

    evidence_supported = bool(
        row.get(
            "biological_evidence_supported"
        )
    )

    evidence_recall = _number_or_none(
        row.get(
            "biological_evidence_recall"
        )
    )

    evidence_jaccard = _number_or_none(
        row.get(
            "biological_evidence_jaccard"
        )
    )

    evidence_precision = _number_or_none(
        row.get(
            "biological_evidence_precision"
        )
    )

    evidence_overlap = _number_or_none(
        row.get(
            "biological_evidence_overlap_count"
        )
    )

    return (
        _pose_validity_priority(row),
        (
            1.0
            if (
                evidence_used
                and evidence_supported
            )
            else 0.0
        ),
        (
            evidence_recall
            if (
                evidence_used
                and evidence_recall
                is not None
            )
            else 0.0
        ),
        (
            evidence_jaccard
            if (
                evidence_used
                and evidence_jaccard
                is not None
            )
            else 0.0
        ),
        (
            evidence_precision
            if (
                evidence_used
                and evidence_precision
                is not None
            )
            else 0.0
        ),
        (
            evidence_overlap
            if (
                evidence_used
                and evidence_overlap
                is not None
            )
            else 0.0
        ),
        (
            cnn_score
            if cnn_score is not None
            else float("-inf")
        ),
        (
            cnn_affinity
            if cnn_affinity is not None
            else float("-inf")
        ),
        (
            -minimized_affinity
            if minimized_affinity is not None
            else float("-inf")
        ),
        (
            -pocket_rank
            if pocket_rank is not None
            else float("-inf")
        ),
    )


def rank_pocket_attempts(
    rows: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    ranked = [
        dict(row)
        for row in rows
    ]

    ranked.sort(
        key=_row_selection_key,
        reverse=True,
    )

    for selection_rank, row in enumerate(
        ranked,
        start=1,
    ):
        row["selection_rank"] = (
            selection_rank
        )
        row["selected"] = (
            selection_rank == 1
        )

    return ranked


def write_pocket_selection_summary(
    output_dir: Path,
    rows: Iterable[dict[str, Any]],
) -> tuple[Path, Path] | None:
    row_list = list(rows)

    if not row_list:
        return None

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    csv_path = (
        output_dir
        / "pocket_selection_summary.csv"
    )
    json_path = (
        output_dir
        / "pocket_selection_summary.json"
    )

    fieldnames = [
        "compound",
        "receptor_conformer_id",
        "selection_rank",
        "selected",
        "pocket_id",
        "pocket_rank",
        "fpocket_score",
        "pocket_source",
        "raw_poses",
        "accepted_poses",
        "rejected_poses",
        "score_source",
        "biological_evidence_available",
        "biological_evidence_used_for_selection",
        "biological_evidence_supported",
        "biological_evidence_rank",
        "biological_evidence_id",
        "biological_evidence_origin",
        "biological_evidence_confidence",
        "biological_evidence_overlap_count",
        "biological_evidence_recall",
        "biological_evidence_precision",
        "biological_evidence_jaccard",
        "biological_evidence_matched_residues",
        "pocket_lining_residue_count",
        "pocket_lining_source_status",
        "top_cnn_score",
        "top_cnn_affinity",
        "top_minimized_affinity",
        "top_seed",
        "top_pose_number",
    ]

    with csv_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()

        for row in row_list:
            serialized = {
                field: row.get(field)
                for field in fieldnames
            }

            for field in (
                "biological_evidence_matched_residues",
            ):
                value = serialized.get(
                    field
                )

                if isinstance(value, list):
                    serialized[field] = (
                        "|".join(value)
                    )

            writer.writerow(serialized)

    selected_rows = [
        row
        for row in row_list
        if row.get("selected") is True
    ]

    evidence_available = any(
        bool(
            row.get(
                "biological_evidence_available"
            )
        )
        for row in row_list
    )

    evidence_used = any(
        bool(
            row.get(
                "biological_evidence_used_for_selection"
            )
        )
        for row in row_list
    )

    if evidence_used:
        selection_method = (
            "PoseBusters-accepted poses preferred; "
            "biologically supported pockets prioritized "
            "by residue-evidence recall, Jaccard, "
            "precision, and overlap; GNINA CNNscore, "
            "CNNaffinity, and minimized affinity used "
            "after biological support"
        )
    else:
        selection_method = (
            "PoseBusters-accepted poses preferred over "
            "raw fallback; highest top-pose CNNscore; "
            "CNNaffinity then minimized affinity used "
            "as tie-breakers"
        )

    first_evidence_row = next(
        (
            row
            for row in row_list
            if row.get(
                "biological_evidence_available"
            )
        ),
        {},
    )

    payload = {
        "selection_method": (
            selection_method
        ),
        "accepted_poses_preferred": True,
        "raw_pose_fallback_when_no_posebusters_pose": True,
        "raw_pose_fallback_ranked_below_accepted_poses": True,
        "biological_evidence_available": (
            evidence_available
        ),
        "biological_evidence_used_for_selection": (
            evidence_used
        ),
        "biological_evidence_id": (
            first_evidence_row.get(
                "biological_evidence_id"
            )
        ),
        "biological_evidence_origin": (
            first_evidence_row.get(
                "biological_evidence_origin"
            )
        ),
        "reference_ligand_used_for_selection": False,
        "compound_count": len(
            {
                str(row.get("compound"))
                for row in row_list
            }
        ),
        "attempt_count": len(row_list),
        "receptor_conformer_count": len(
            {
                str(
                    row.get(
                        "receptor_conformer_id",
                        "submitted_receptor",
                    )
                )
                for row in row_list
            }
        ),
        "selection_unit": (
            "receptor_conformer_pocket_pair"
        ),
        "selected_pockets": selected_rows,
        "attempts": row_list,
    }

    json_path.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return csv_path, json_path
