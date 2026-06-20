from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable

from .models import PocketDefinition, PoseRecord


def _number_or_none(value: object) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)

    return None


def _pose_selection_key(
    record: PoseRecord,
) -> tuple[float, float, float]:
    cnn_score = _number_or_none(record.cnn_score)
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
    pocket: PocketDefinition,
    raw_records: Iterable[PoseRecord],
    accepted_records: Iterable[PoseRecord],
    rejected_pose_count: int,
) -> dict[str, Any]:
    raw_list = list(raw_records)
    accepted_list = list(accepted_records)

    if accepted_list:
        scoring_records = accepted_list
        score_source = "accepted_poses"
    elif raw_list:
        scoring_records = raw_list
        score_source = "raw_pose_fallback"
    else:
        scoring_records = []
        score_source = "none"

    best_record = best_pose_record(scoring_records)

    return {
        "compound": ligand_name,
        "pocket_id": pocket.pocket_id,
        "pocket_rank": pocket.pocket_rank,
        "fpocket_score": pocket.fpocket_score,
        "pocket_source": (
            pocket.source or pocket.mode
        ),
        "raw_poses": len(raw_list),
        "accepted_poses": len(accepted_list),
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


def _row_selection_key(
    row: dict[str, Any],
) -> tuple[float, float, float, float]:
    cnn_score = _number_or_none(
        row.get("top_cnn_score")
    )
    cnn_affinity = _number_or_none(
        row.get("top_cnn_affinity")
    )
    minimized_affinity = _number_or_none(
        row.get("top_minimized_affinity")
    )
    pocket_rank = _number_or_none(
        row.get("pocket_rank")
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
        row["selection_rank"] = selection_rank
        row["selected"] = selection_rank == 1

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
            writer.writerow(
                {
                    field: row.get(field)
                    for field in fieldnames
                }
            )

    selected_rows = [
        row
        for row in row_list
        if row.get("selected") is True
    ]

    payload = {
        "selection_method": (
            "highest top-pose CNNscore; "
            "CNNaffinity then minimized affinity "
            "used as tie-breakers"
        ),
        "accepted_poses_preferred": True,
        "raw_pose_fallback_when_no_posebusters_pose": True,
        "reference_ligand_used_for_selection": False,
        "compound_count": len(
            {
                str(row.get("compound"))
                for row in row_list
            }
        ),
        "attempt_count": len(row_list),
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
