from __future__ import annotations

from .models import PoseCluster


def _score(cluster: PoseCluster) -> float | None:
    return cluster.representative.cnn_score


def assess_uncertainty(
    clusters: list[PoseCluster],
    seed_count: int,
) -> tuple[str, list[str]]:
    if not clusters:
        return "unresolved", ["No physically valid pose hypotheses remained."]

    reasons: list[str] = []
    top = clusters[0]
    top_score = _score(top)
    second_score = _score(clusters[1]) if len(clusters) > 1 else None
    top_seed_fraction = len(top.seeds) / max(seed_count, 1)
    top_population = top.member_count / max(sum(c.member_count for c in clusters), 1)

    if second_score is None or top_score is None or second_score is None:
        score_gap = None
    else:
        score_gap = top_score - second_score

    if len(clusters) == 1:
        reasons.append("Only one physically valid pose cluster remained.")
    else:
        reasons.append(f"{len(clusters)} distinct physically valid pose clusters remained.")

    reasons.append(
        f"The leading cluster appeared in {len(top.seeds)}/{max(seed_count, 1)} seeds "
        f"and contained {top_population:.0%} of valid poses."
    )
    if score_gap is not None:
        reasons.append(f"GNINA CNN score gap to the second cluster was {score_gap:.3f}.")

    if (
        len(clusters) == 1
        or (
            top_seed_fraction >= 0.67
            and top_population >= 0.50
            and score_gap is not None
            and score_gap >= 0.15
        )
    ):
        return "high", reasons

    if (
        top_seed_fraction >= 0.50
        and score_gap is not None
        and score_gap >= 0.08
    ):
        return "moderate", reasons

    return "low", reasons
