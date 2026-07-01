"""Production-style evidence-aware pose retention from final hypothesis PDB files.

This module does not use RMSD. It reads final hypothesis PDB REMARK 900 metadata
and scores alternatives by:

- GNINA CNN score
- contact overlap with optional known/expected residues
- contact recurrence across hypotheses
- cluster/seed support
- physical sanity from closest contact and minimized affinity

Policy:
- Keep the top-CNN pose as the primary hypothesis.
- Retain alternatives when they show stronger contact, recurrence, clustering,
  or physical-sanity evidence.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ContactRetentionConfig:
    near_cnn_delta: float = 0.10
    min_recurrent_contacts: int = 2
    max_retained_per_compound: int = 5


def parse_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null", "unknown", "unavailable"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def normalize_residue(token: str) -> str:
    token = token.strip().strip(",;")
    token = token.replace(" ", "")
    return token.upper()


def parse_residue_list(text: str) -> set[str]:
    if not text:
        return set()

    # Match common forms like ARG118:B, ASP151:B, GLY147:B.
    matches = re.findall(r"\b[A-Z]{3}\d+[A-Z]?:[A-Za-z0-9]+\b", text.upper())
    return {normalize_residue(item) for item in matches}


def split_known_residues(raw: str | None) -> set[str]:
    if not raw:
        return set()
    parts = re.split(r"[,\s;]+", raw)
    return {normalize_residue(part) for part in parts if part.strip()}


def parse_remark_metadata(pdb_path: Path) -> dict:
    metadata = {
        "pdb_file": str(pdb_path),
        "pdb_name": pdb_path.name,
        "compound": "",
        "compound_priority_rank": "",
        "hypothesis_index": "",
        "hypothesis_total": "",
        "cnn_score": None,
        "cnn_affinity": None,
        "minimized_affinity": None,
        "cluster_members": None,
        "seeds_represented": None,
        "pose_confidence": "",
        "pocket_id": "",
        "pocket_source": "",
        "receptor_conformer": "",
        "physical_validity": "",
        "closest_contact_angstrom": None,
        "contact_residues": set(),
        "polar_contacts": set(),
        "hydrophobic_contacts": set(),
    }

    last_key = None

    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith("REMARK 900"):
            continue

        payload = line.replace("REMARK 900", "", 1).strip()

        if payload.startswith("COMPOUND "):
            metadata["compound"] = payload.replace("COMPOUND ", "", 1).strip()
            last_key = "compound"

        elif payload.startswith("COMPOUND PRIORITY RANK"):
            metadata["compound_priority_rank"] = payload.replace(
                "COMPOUND PRIORITY RANK", "", 1
            ).strip()
            last_key = None

        elif payload.startswith("HYPOTHESIS"):
            match = re.search(r"HYPOTHESIS\s+(\d+)\s+OF\s+(\d+)", payload)
            if match:
                metadata["hypothesis_index"] = match.group(1)
                metadata["hypothesis_total"] = match.group(2)
            last_key = None

        elif payload.startswith("GNINA CNN SCORE"):
            metadata["cnn_score"] = parse_float(
                payload.replace("GNINA CNN SCORE", "", 1)
            )
            last_key = None

        elif payload.startswith("GNINA CNN AFFINITY"):
            metadata["cnn_affinity"] = parse_float(
                payload.replace("GNINA CNN AFFINITY", "", 1)
            )
            last_key = None

        elif payload.startswith("GNINA MINIMIZED AFFINITY"):
            metadata["minimized_affinity"] = parse_float(
                payload.replace("GNINA MINIMIZED AFFINITY", "", 1)
            )
            last_key = None

        elif payload.startswith("CLUSTER MEMBERS"):
            metadata["cluster_members"] = parse_float(
                payload.replace("CLUSTER MEMBERS", "", 1)
            )
            last_key = None

        elif payload.startswith("SEEDS REPRESENTED"):
            metadata["seeds_represented"] = parse_float(
                payload.replace("SEEDS REPRESENTED", "", 1)
            )
            last_key = None

        elif payload.startswith("POSE CONFIDENCE"):
            metadata["pose_confidence"] = payload.replace("POSE CONFIDENCE", "", 1).strip()
            last_key = None

        elif payload.startswith("POCKET ID"):
            metadata["pocket_id"] = payload.replace("POCKET ID", "", 1).strip()
            last_key = None

        elif payload.startswith("POCKET SOURCE"):
            metadata["pocket_source"] = payload.replace("POCKET SOURCE", "", 1).strip()
            last_key = None

        elif payload.startswith("RECEPTOR CONFORMER"):
            metadata["receptor_conformer"] = payload.replace(
                "RECEPTOR CONFORMER", "", 1
            ).strip()
            last_key = None

        elif payload.startswith("PHYSICAL VALIDITY"):
            metadata["physical_validity"] = payload.replace(
                "PHYSICAL VALIDITY", "", 1
            ).strip()
            last_key = None

        elif payload.startswith("CLOSEST PROTEIN CONTACT"):
            match = re.search(r"([0-9.]+)\s+ANGSTROM", payload)
            if match:
                metadata["closest_contact_angstrom"] = parse_float(match.group(1))
            last_key = None

        elif payload.startswith("CONTACT RESIDUES"):
            residues = parse_residue_list(payload.replace("CONTACT RESIDUES", "", 1))
            metadata["contact_residues"].update(residues)
            last_key = "contact_residues"

        elif payload.startswith("POLAR CONTACT CANDIDATES"):
            residues = parse_residue_list(
                payload.replace("POLAR CONTACT CANDIDATES", "", 1)
            )
            metadata["polar_contacts"].update(residues)
            last_key = "polar_contacts"

        elif payload.startswith("HYDROPHOBIC CONTACT RESIDUES"):
            residues = parse_residue_list(
                payload.replace("HYDROPHOBIC CONTACT RESIDUES", "", 1)
            )
            metadata["hydrophobic_contacts"].update(residues)
            last_key = "hydrophobic_contacts"

        elif last_key in {"contact_residues", "polar_contacts", "hydrophobic_contacts"}:
            metadata[last_key].update(parse_residue_list(payload))

    if str(metadata["compound"]).upper().startswith("PRIORITY RANK"):
        metadata["compound"] = ""

    if not metadata["compound"]:
        # Fallback from filename: 01__oseltamivir__hypothesis_01.pdb
        parts = pdb_path.stem.split("__")
        if len(parts) >= 2:
            metadata["compound"] = parts[1]
        else:
            metadata["compound"] = pdb_path.stem

    return metadata


def physical_sanity_score(pose: dict) -> tuple[float, list[str]]:
    score = 0.0
    tags = []

    closest = pose.get("closest_contact_angstrom")
    if closest is not None:
        if 1.5 <= closest <= 4.0:
            score += 1.0
            tags.append("sane_contact_distance")
        elif closest < 1.2:
            score -= 1.0
            tags.append("possible_clash")
        elif closest > 5.0:
            score -= 0.5
            tags.append("weak_contact_distance")

    minimized_affinity = pose.get("minimized_affinity")
    if minimized_affinity is not None:
        if minimized_affinity < 50:
            score += 0.5
            tags.append("sane_minimized_affinity")
        else:
            score -= 1.0
            tags.append("extreme_minimized_affinity")

    validity = str(pose.get("physical_validity") or "").upper()
    if "PASS" in validity:
        score += 1.0
        tags.append("physical_validity_pass")

    return score, tags


def score_pose(
    pose: dict,
    *,
    primary_pose: dict | None,
    recurrent_residues: set[str],
    known_residues: set[str],
) -> dict:
    contacts = set(pose.get("contact_residues") or set())

    cnn = pose.get("cnn_score")
    cnn_component = cnn if cnn is not None else 0.0

    known_overlap = contacts & known_residues if known_residues else set()
    recurrent_overlap = contacts & recurrent_residues if recurrent_residues else set()

    primary_overlap = set()
    if primary_pose:
        primary_contacts = set(primary_pose.get("contact_residues") or set())
        primary_overlap = contacts & primary_contacts

    sanity_score, sanity_tags = physical_sanity_score(pose)

    cluster_members = pose.get("cluster_members") or 0.0
    seeds_represented = pose.get("seeds_represented") or 0.0

    evidence_score = (
        cnn_component
        + 3.0 * len(known_overlap)
        + 1.0 * len(recurrent_overlap)
        + 0.25 * len(primary_overlap)
        + 0.25 * cluster_members
        + 0.50 * seeds_represented
        + sanity_score
    )

    tags = []
    if known_overlap:
        tags.append("known_residue_overlap")
    if recurrent_overlap:
        tags.append("recurrent_contact_overlap")
    if primary_overlap:
        tags.append("primary_contact_overlap")
    tags.extend(sanity_tags)

    return {
        "evidence_score": evidence_score,
        "known_overlap": known_overlap,
        "recurrent_overlap": recurrent_overlap,
        "primary_overlap": primary_overlap,
        "physical_sanity_score": sanity_score,
        "evidence_tags": tags,
    }


def analyze_output_dir(
    output_dir: Path,
    *,
    known_residues: set[str] | None = None,
    config: ContactRetentionConfig | None = None,
):
    output_dir = Path(output_dir)
    known_residues = known_residues or set()
    config = config or ContactRetentionConfig()

    hypothesis_files = sorted(output_dir.glob("*hypothesis*.pdb"))

    if not hypothesis_files:
        raise FileNotFoundError(f"No *hypothesis*.pdb files found in {output_dir}")

    poses = [parse_remark_metadata(path) for path in hypothesis_files]

    by_compound: dict[str, list[dict]] = defaultdict(list)
    for pose in poses:
        by_compound[pose["compound"]].append(pose)

    rows = []

    for compound, compound_poses in sorted(by_compound.items()):
        compound_poses = sorted(
            compound_poses,
            key=lambda pose: (
                -(pose["cnn_score"] if pose["cnn_score"] is not None else -1e9),
                pose["pdb_name"],
            ),
        )

        for rank, pose in enumerate(compound_poses, start=1):
            pose["cnn_rank_within_compound"] = rank

        primary_pose = compound_poses[0]

        contact_counter = Counter()
        for pose in compound_poses:
            contact_counter.update(set(pose.get("contact_residues") or set()))

        recurrent_residues = {
            residue
            for residue, count in contact_counter.items()
            if count >= config.min_recurrent_contacts
        }

        scored = []
        for pose in compound_poses:
            evidence = score_pose(
                pose,
                primary_pose=primary_pose,
                recurrent_residues=recurrent_residues,
                known_residues=known_residues,
            )
            scored.append((pose, evidence))

        primary_evidence = next(
            evidence for pose, evidence in scored if pose is primary_pose
        )

        top_cnn = primary_pose.get("cnn_score")
        primary_evidence_score = primary_evidence["evidence_score"]

        for pose, evidence in scored:
            is_primary = pose is primary_pose

            cnn_delta = ""
            if top_cnn is not None and pose.get("cnn_score") is not None:
                cnn_delta = top_cnn - pose["cnn_score"]

            retain_reasons = []
            candidate_type = "supporting_alternative"

            if is_primary:
                retain_reasons.append("primary_top_cnn")
                candidate_type = "primary_top_cnn"

            if evidence["known_overlap"]:
                retain_reasons.append("known_contact_supported")
            if evidence["recurrent_overlap"]:
                retain_reasons.append("ensemble_recurrent_contacts")
            if evidence["evidence_score"] > primary_evidence_score and not is_primary:
                retain_reasons.append("evidence_score_exceeds_primary")
            if (
                cnn_delta != ""
                and cnn_delta <= config.near_cnn_delta
                and not is_primary
            ):
                retain_reasons.append("near_top_cnn")

            retain = bool(retain_reasons)

            rows.append(
                {
                    "compound": compound,
                    "candidate_type": candidate_type,
                    "retain_recommended": retain,
                    "hypothesis_file": pose["pdb_name"],
                    "cnn_rank_within_compound": pose["cnn_rank_within_compound"],
                    "cnn_score": pose.get("cnn_score", ""),
                    "cnn_delta_from_top": cnn_delta,
                    "evidence_score": evidence["evidence_score"],
                    "primary_evidence_score": primary_evidence_score,
                    "physical_sanity_score": evidence["physical_sanity_score"],
                    "cluster_members": pose.get("cluster_members", ""),
                    "seeds_represented": pose.get("seeds_represented", ""),
                    "closest_contact_angstrom": pose.get("closest_contact_angstrom", ""),
                    "minimized_affinity": pose.get("minimized_affinity", ""),
                    "known_contact_overlap_count": len(evidence["known_overlap"]),
                    "known_contact_overlap": ";".join(sorted(evidence["known_overlap"])),
                    "recurrent_contact_overlap_count": len(evidence["recurrent_overlap"]),
                    "recurrent_contact_overlap": ";".join(sorted(evidence["recurrent_overlap"])),
                    "primary_contact_overlap_count": len(evidence["primary_overlap"]),
                    "primary_contact_overlap": ";".join(sorted(evidence["primary_overlap"])),
                    "all_contacts": ";".join(sorted(pose.get("contact_residues") or set())),
                    "evidence_tags": ";".join(evidence["evidence_tags"]),
                    "retain_reasons": ";".join(retain_reasons),
                    "interpretation": interpret_contact_retention_row(
                        is_primary=is_primary,
                        retain_reasons=retain_reasons,
                        evidence=evidence,
                    ),
                }
            )

    return rows


def interpret_contact_retention_row(*, is_primary: bool, retain_reasons: list[str], evidence: dict) -> str:
    if is_primary:
        return "Primary GNINA/CNN hypothesis retained as the main pose."

    if "evidence_score_exceeds_primary" in retain_reasons:
        return "Alternative has stronger total evidence score than the top-CNN pose."

    if "known_contact_supported" in retain_reasons:
        return "Alternative overlaps known/expected contact residues."

    if "ensemble_recurrent_contacts" in retain_reasons:
        return "Alternative shares recurrent contact residues seen across hypotheses."

    if "near_top_cnn" in retain_reasons:
        return "Alternative is close to the top CNN score and retained for review."

    return "Alternative not retained by current evidence policy."


def write_contact_retention_artifacts(
    output_dir: Path,
    *,
    known_residues: set[str] | None = None,
    config: ContactRetentionConfig | None = None,
):
    output_dir = Path(output_dir)
    known_residues = known_residues or set()
    config = config or ContactRetentionConfig()

    rows = analyze_output_dir(
        output_dir,
        known_residues=known_residues,
        config=config,
    )

    out_csv = output_dir / "production_pose_retention_candidates.csv"
    out_md = output_dir / "production_pose_retention_report.md"

    fieldnames = [
        "compound",
        "candidate_type",
        "retain_recommended",
        "hypothesis_file",
        "cnn_rank_within_compound",
        "cnn_score",
        "cnn_delta_from_top",
        "evidence_score",
        "primary_evidence_score",
        "physical_sanity_score",
        "cluster_members",
        "seeds_represented",
        "closest_contact_angstrom",
        "minimized_affinity",
        "known_contact_overlap_count",
        "known_contact_overlap",
        "recurrent_contact_overlap_count",
        "recurrent_contact_overlap",
        "primary_contact_overlap_count",
        "primary_contact_overlap",
        "all_contacts",
        "evidence_tags",
        "retain_reasons",
        "interpretation",
    ]

    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    retained = [row for row in rows if str(row["retain_recommended"]) == "True"]

    md = [
        "# Production Pose Retention Report",
        "",
        f"- Output directory: `{output_dir}`",
        f"- Hypothesis poses analyzed: {len(rows)}",
        f"- Retained poses: {len(retained)}",
        f"- Optional known residues supplied: {', '.join(sorted(known_residues)) if known_residues else 'none'}",
        "",
        "## Policy",
        "",
        "The top GNINA CNN pose remains the primary hypothesis. Alternative poses are retained when they show stronger production-style evidence such as contact overlap, contact recurrence, physical sanity, clustering support, or near-top CNN score.",
        "",
        "## Retained / Scored Poses",
        "",
        "| Compound | Retain | Type | File | CNN rank | CNN score | Evidence score | Known overlap | Recurrent overlap | Reasons | Interpretation |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---|---|",
    ]

    for row in rows:
        md.append(
            "| {compound} | {retain_recommended} | {candidate_type} | `{hypothesis_file}` | {cnn_rank_within_compound} | {cnn_score} | {evidence_score:.3f} | {known_contact_overlap_count} | {recurrent_contact_overlap_count} | {retain_reasons} | {interpretation} |".format(
                **row
            )
        )

    md.extend(
        [
            "",
            "## Notes",
            "",
            "- This report does not use RMSD and is suitable for production-style runs.",
            "- Contact recurrence is computed within the available hypotheses in this output directory.",
            "- Known residue overlap is only used when expected/active-site residues are supplied.",
            "- CNN score is retained as one signal, not the only signal.",
            "",
            f"CSV output: `{out_csv.name}`",
            "",
        ]
    )

    out_md.write_text("\n".join(md), encoding="utf-8")

    return out_csv, out_md


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--known-residues",
        default="",
        help="Comma/space separated expected residues, e.g. ARG118:B,ASP151:B",
    )
    parser.add_argument("--near-cnn-delta", type=float, default=0.10)
    parser.add_argument("--min-recurrent-contacts", type=int, default=2)
    parser.add_argument("--max-retained-per-compound", type=int, default=5)
    args = parser.parse_args(argv)

    config = ContactRetentionConfig(
        near_cnn_delta=args.near_cnn_delta,
        min_recurrent_contacts=args.min_recurrent_contacts,
        max_retained_per_compound=args.max_retained_per_compound,
    )

    known_residues = split_known_residues(args.known_residues)

    out_csv, out_md = write_contact_retention_artifacts(
        args.output_dir,
        known_residues=known_residues,
        config=config,
    )

    print("Wrote:", out_csv)
    print("Wrote:", out_md)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
