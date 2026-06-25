from __future__ import annotations

from pathlib import Path

from .models import PoseRecord


def receptor_display_pdb_for_pose(
    record: PoseRecord,
    fallback_receptor_pdb: Path,
) -> Path:
    """Return the receptor coordinates matching a pose."""

    if record.receptor_display_pdb is not None:
        return Path(
            record.receptor_display_pdb
        )

    return Path(
        fallback_receptor_pdb
    )


def records_share_receptor_conformer(
    first: PoseRecord,
    second: PoseRecord,
) -> bool:
    """Whether two poses came from the same receptor."""

    return (
        first.receptor_conformer_id
        == second.receptor_conformer_id
    )
