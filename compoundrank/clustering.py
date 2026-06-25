from __future__ import annotations

from .chemistry import direct_heavy_rmsd
from .conformer_context import (
    records_share_receptor_conformer,
)
from .models import PoseCluster, PoseRecord


def gnina_sort_key(record: PoseRecord) -> tuple[float, float, int, int]:
    cnn_score = record.cnn_score if record.cnn_score is not None else float("-inf")
    cnn_affinity = (
        record.cnn_affinity if record.cnn_affinity is not None else float("-inf")
    )
    return (cnn_score, cnn_affinity, -record.seed, -record.pose_number)


def cluster_pose_hypotheses(
    records: list[PoseRecord],
    *,
    rmsd_threshold: float = 2.0,
) -> list[PoseCluster]:
    """Cluster poses without inventing a new score.

    Poses are processed in GNINA score order. The first pose that establishes a
    cluster remains its representative. Cluster ordering therefore remains the
    ordering provided by GNINA, while geometrically redundant poses are grouped.
    """
    ordered = sorted(records, key=gnina_sort_key, reverse=True)
    clusters: list[PoseCluster] = []

    for record in ordered:
        assigned_cluster: PoseCluster | None = None
        for cluster in clusters:
            if not records_share_receptor_conformer(
                record,
                cluster.representative,
            ):
                continue

            rmsd = direct_heavy_rmsd(
                record.molecule,
                cluster.representative.molecule,
            )

            if rmsd <= rmsd_threshold:
                assigned_cluster = cluster
                break

        if assigned_cluster is None:
            assigned_cluster = PoseCluster(
                cluster_id=len(clusters) + 1,
                representative=record,
                members=[],
                seeds=set(),
                member_count=0,
                valid_member_count=0,
            )
            clusters.append(assigned_cluster)

        assigned_cluster.members.append(record)
        assigned_cluster.seeds.add(record.seed)
        assigned_cluster.member_count += 1
        assigned_cluster.valid_member_count += 1

    return clusters
