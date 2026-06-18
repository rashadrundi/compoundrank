"""Stage 4A: evidence-grounded compound retrieval.

This module turns target evidence into candidate ligand packets.

Current stage:
- local rule-based retrieval
- optional PubChem structure fetch by compound name
- candidate JSON/CSV output
- docking manifest output
- ligand search report output

Future expansion hooks:
- ChEMBL / BindingDB / RCSB retrieval
- mutation-aware ligand adaptation
- pocket-feature interpretation
- novel ligand hypothesis generation
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from .ligand_rules import LIGAND_RETRIEVAL_RULES


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def flatten_text(value: Any) -> str:
    """Flatten nested JSON-like data into searchable lowercase text."""
    parts: list[str] = []

    def visit(item: Any) -> None:
        if item is None:
            return
        if isinstance(item, dict):
            for key, val in item.items():
                parts.append(str(key))
                visit(val)
            return
        if isinstance(item, (list, tuple, set)):
            for element in item:
                visit(element)
            return
        parts.append(str(item))

    visit(value)
    return " ".join(parts).lower()


def _nested(value: dict[str, Any], *keys: str) -> Any:
    current: Any = value
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def _first_text(*values: Any, default: str = "unknown") -> str:
    for value in values:
        if value is None:
            continue
        if isinstance(value, (dict, list, tuple, set)):
            continue
        text = str(value).strip()
        if text:
            return text
    return default


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value if str(item).strip()]
    if isinstance(value, tuple):
        return [str(item) for item in value if str(item).strip()]
    if isinstance(value, str) and value.strip():
        return [value.strip()]
    return []


def target_context(target_evidence: dict[str, Any]) -> dict[str, Any]:
    """Normalize target-evidence fields across output schema versions."""
    interpretation = target_evidence.get("target_interpretation") or {}
    special_domain = _nested(target_evidence, "evidence", "special_domain_evidence") or {}

    query_terms = (
        target_evidence.get("future_ligand_database_query")
        or target_evidence.get("future_ligand_database_query_terms")
        or target_evidence.get("ligand_query_terms")
        or interpretation.get("future_ligand_database_query")
        or interpretation.get("future_ligand_database_query_terms")
        or interpretation.get("ligand_query_terms")
        or []
    )

    return {
        "target_name": _first_text(
            target_evidence.get("target_name"),
            interpretation.get("target_name"),
            special_domain.get("target_name"),
        ),
        "target_class": _first_text(
            target_evidence.get("target_class"),
            interpretation.get("target_class"),
            special_domain.get("target_class"),
        ),
        "enzyme_class": _first_text(
            target_evidence.get("enzyme_class"),
            interpretation.get("enzyme_class"),
            special_domain.get("enzyme_class"),
        ),
        "viral_family": _first_text(
            target_evidence.get("viral_family_evidence"),
            interpretation.get("viral_family_evidence"),
            interpretation.get("viral_family"),
            special_domain.get("viral_family"),
        ),
        "confidence": _first_text(
            target_evidence.get("evidence_confidence"),
            target_evidence.get("confidence"),
            interpretation.get("evidence_confidence"),
            interpretation.get("confidence"),
            special_domain.get("confidence"),
        ),
        "special_domain_label": _first_text(
            special_domain.get("label"),
            special_domain.get("hit_name"),
            default="",
        ),
        "special_domain_accession": _first_text(
            special_domain.get("accession"),
            default="",
        ),
        "query_terms": _as_list(query_terms),
    }


def matching_terms(rule: dict[str, Any], evidence_blob: str) -> list[str]:
    matches: list[str] = []
    for term in rule.get("match_terms", []):
        if str(term).lower() in evidence_blob:
            matches.append(str(term))
    return matches


def candidate_key(candidate: dict[str, Any]) -> str:
    return str(candidate.get("compound_name", "")).strip().lower()


def make_candidate(
    *,
    compound: dict[str, Any],
    rule: dict[str, Any],
    matched_terms: list[str],
    target_evidence: dict[str, Any],
    retrieval_rank: int,
) -> dict[str, Any]:
    compound_name = str(compound["name"])
    target = target_context(target_evidence)

    evidence_level = str(
        compound.get("evidence_level")
        or rule.get("rule_evidence_level")
        or "moderate"
    )

    retrieval_reason = (
        f"Selected by rule '{rule.get('rule_id')}' because target evidence matched "
        f"{', '.join(matched_terms) if matched_terms else 'rule terms'}."
    )

    candidate = {
        "compound_name": compound_name,
        "retrieval_rank": retrieval_rank,
        "design_status": compound.get("design_status", "known_inhibitor"),
        "retrieval_route": "domain_family_rule",
        "retrieval_rule_id": rule.get("rule_id"),
        "retrieval_rule_label": rule.get("rule_label"),
        "retrieval_reason": retrieval_reason,
        "target_family_basis": rule.get("target_family_basis"),
        "domain_basis": matched_terms,
        "target_name": target["target_name"],
        "target_class": target["target_class"] if target["target_class"] != "unknown" else rule.get("target_class"),
        "enzyme_class": target["enzyme_class"] if target["enzyme_class"] != "unknown" else rule.get("enzyme_class"),
        "viral_family": target["viral_family"],
        "target_evidence_confidence": target["confidence"],
        "special_domain_label": target["special_domain_label"],
        "special_domain_accession": target["special_domain_accession"],
        "evidence_level": evidence_level,
        "retrieval_terms": rule.get("retrieval_terms", []),
        "source_databases": ["local_rule_registry"],
        "source_notes": compound.get("notes", ""),
        "pubchem_cid": None,
        "smiles": None,
        "inchi_key": None,
        "structure_source": None,
        "structure_fetch_status": "not_attempted",
        "structure_fetch_error": None,
        "local_sdf_path": None,
        "selected_for_docking": True,
        "mutation_relevance": {
            "status": "not_evaluated",
            "pocket_mutation_overlap": None,
            "expected_interaction_change": None,
            "adaptation_rationale": None,
        },
        "novel_design_relevance": {
            "status": "not_generated",
            "pocket_feature_basis": None,
            "pharmacophore_basis": None,
            "hypothetical_modification": None,
        },
        "limitations": rule.get("limitations", []),
    }

    return candidate


def build_ligand_candidates(
    target_evidence: dict[str, Any],
    homolog_summary: dict[str, Any] | None = None,
    *,
    max_candidates: int = 20,
) -> list[dict[str, Any]]:
    """Build candidate ligand packets from target evidence and optional homology summary."""
    evidence_blob = flatten_text(
        {"target_evidence": target_evidence, "homolog_summary": homolog_summary or {}}
    )

    candidates_by_key: dict[str, dict[str, Any]] = {}

    for rule in LIGAND_RETRIEVAL_RULES:
        terms = matching_terms(rule, evidence_blob)
        if not terms:
            continue

        for retrieval_rank, compound in enumerate(rule.get("seed_compounds", []), start=1):
            candidate = make_candidate(
                compound=compound,
                rule=rule,
                matched_terms=terms,
                target_evidence=target_evidence,
                retrieval_rank=retrieval_rank,
            )
            key = candidate_key(candidate)
            if not key:
                continue

            existing = candidates_by_key.get(key)
            if existing is None:
                candidates_by_key[key] = candidate
                continue

            old_level = str(existing.get("evidence_level", "moderate"))
            new_level = str(candidate.get("evidence_level", "moderate"))
            priority = {"strong": 3, "moderate": 2, "weak": 1, "exploratory": 0}
            if priority.get(new_level, 0) > priority.get(old_level, 0):
                candidates_by_key[key] = candidate

    candidates = list(candidates_by_key.values())

    evidence_priority = {"strong": 0, "moderate": 1, "weak": 2, "exploratory": 3}
    candidates.sort(
        key=lambda item: (
            evidence_priority.get(str(item.get("evidence_level", "moderate")), 1),
            int(item.get("retrieval_rank") or 9999),
            str(item.get("compound_name", "")).lower(),
        )
    )

    if max_candidates > 0:
        candidates = candidates[:max_candidates]

    return candidates


def pubchem_request_text(url: str, *, timeout_seconds: int) -> str:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "CompoundRank-stage4A/0.1 educational-research-use"},
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        return response.read().decode("utf-8", errors="replace").strip()


def pubchem_request_bytes(url: str, *, timeout_seconds: int) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "CompoundRank-stage4A/0.1 educational-research-use"},
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        return response.read()


def fetch_pubchem_cid_by_name(name: str, *, timeout_seconds: int) -> str | None:
    encoded = urllib.parse.quote(name)
    url = f"{PUBCHEM_BASE}/compound/name/{encoded}/cids/TXT"
    text = pubchem_request_text(url, timeout_seconds=timeout_seconds)
    first = text.splitlines()[0].strip() if text else ""
    return first or None


def fetch_pubchem_sdf_by_name(
    name: str,
    destination: Path,
    *,
    timeout_seconds: int,
) -> tuple[str, str]:
    """Fetch an SDF from PubChem by name.

    Returns:
        (record_type, local_path)

    Tries 3D first, then 2D.
    """
    encoded = urllib.parse.quote(name)

    errors: list[str] = []
    for record_type in ("3d", "2d"):
        url = f"{PUBCHEM_BASE}/compound/name/{encoded}/SDF?record_type={record_type}"
        try:
            data = pubchem_request_bytes(url, timeout_seconds=timeout_seconds)
            if len(data) < 100:
                raise ValueError("PubChem returned unexpectedly small SDF response.")
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
            return record_type, str(destination)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{record_type}: {exc}")

    raise RuntimeError("; ".join(errors))


def fetch_candidate_structures_from_pubchem(
    candidates: list[dict[str, Any]],
    *,
    output_dir: Path,
    timeout_seconds: int = 60,
    sleep_seconds: float = 0.2,
) -> list[dict[str, Any]]:
    """Fetch candidate SDF structures from PubChem by compound name."""
    ligand_dir = output_dir / "retrieved_ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)

    for candidate in candidates:
        if not candidate.get("selected_for_docking", True):
            continue

        name = str(candidate.get("compound_name", "")).strip()
        if not name:
            continue

        safe_name = (
            name.lower()
            .replace(" ", "_")
            .replace("/", "_")
            .replace("\\", "_")
            .replace(":", "_")
        )
        sdf_path = ligand_dir / f"{safe_name}.sdf"

        try:
            cid = fetch_pubchem_cid_by_name(name, timeout_seconds=timeout_seconds)
            if cid:
                candidate["pubchem_cid"] = cid

            record_type, local_path = fetch_pubchem_sdf_by_name(
                name,
                sdf_path,
                timeout_seconds=timeout_seconds,
            )
            candidate["local_sdf_path"] = local_path
            candidate["structure_source"] = f"PubChem name lookup ({record_type.upper()} SDF)"
            candidate["structure_fetch_status"] = "fetched"
            candidate["structure_fetch_error"] = None

            if "PubChem" not in candidate["source_databases"]:
                candidate["source_databases"].append("PubChem")

        except Exception as exc:  # noqa: BLE001
            candidate["structure_fetch_status"] = "failed"
            candidate["structure_fetch_error"] = str(exc)
            candidate["local_sdf_path"] = None

        time.sleep(sleep_seconds)

    return candidates


def write_candidate_csv(path: Path, candidates: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "compound_name",
        "retrieval_rank",
        "design_status",
        "evidence_level",
        "retrieval_route",
        "retrieval_rule_id",
        "target_family_basis",
        "target_name",
        "target_class",
        "enzyme_class",
        "viral_family",
        "special_domain_label",
        "special_domain_accession",
        "pubchem_cid",
        "structure_fetch_status",
        "local_sdf_path",
        "selected_for_docking",
        "retrieval_reason",
    ]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for candidate in candidates:
            writer.writerow({field: candidate.get(field) for field in fields})


def write_docking_manifest(path: Path, candidates: list[dict[str, Any]]) -> None:
    """Write a ligand manifest compatible with the main CompoundRank pipeline.

    The existing ligand loader requires:
        name, source_type, value

    For PubChem-fetched Stage 4A ligands:
        source_type = file
        value = local SDF path

    Extra metadata columns are preserved for traceability.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "name",
        "source_type",
        "value",
        "retrieval_reason",
        "evidence_level",
        "design_status",
        "pubchem_cid",
    ]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for candidate in candidates:
            local_sdf_path = candidate.get("local_sdf_path")
            if not local_sdf_path:
                continue
            if not candidate.get("selected_for_docking", True):
                continue

            writer.writerow(
                {
                    "name": candidate.get("compound_name"),
                    "source_type": "file",
                    "value": local_sdf_path,
                    "retrieval_reason": candidate.get("retrieval_reason"),
                    "evidence_level": candidate.get("evidence_level"),
                    "design_status": candidate.get("design_status"),
                    "pubchem_cid": candidate.get("pubchem_cid"),
                }
            )


def render_ligand_search_report(
    *,
    target_evidence: dict[str, Any],
    candidates: list[dict[str, Any]],
    fetch_structures: bool,
) -> str:
    lines: list[str] = []
    target = target_context(target_evidence)

    lines.append("# Ligand Search Report")
    lines.append("")
    lines.append("## Target Basis")
    lines.append("")
    lines.append(f"- Target name: {target['target_name']}")
    lines.append(f"- Target class: {target['target_class']}")
    lines.append(f"- Enzyme class: {target['enzyme_class']}")
    lines.append(f"- Viral family evidence: {target['viral_family']}")
    lines.append(f"- Evidence confidence: {target['confidence']}")

    if target["special_domain_label"] or target["special_domain_accession"]:
        lines.append(
            f"- Special domain evidence: {target['special_domain_label']} "
            f"({target['special_domain_accession']})".strip()
        )

    if target["query_terms"]:
        lines.append("- Existing target-evidence query terms:")
        for term in target["query_terms"]:
            lines.append(f"  - {term}")

    lines.append("")
    lines.append("## Retrieval Summary")
    lines.append("")
    lines.append(f"- Candidate count: {len(candidates)}")
    lines.append(f"- Structure fetch attempted: {'yes' if fetch_structures else 'no'}")
    lines.append("")

    if not candidates:
        lines.append("No ligand retrieval rule matched this target evidence.")
        lines.append("")
        lines.append(
            "Future expansion should use pocket-feature interpretation, broader target-family rules, "
            "or novel-ligand hypothesis generation."
        )
        lines.append("")
    else:
        lines.append("| Rank | Compound | Status | Evidence | Rule | Structure |")
        lines.append("|---:|---|---|---|---|---|")
        for candidate in candidates:
            structure_status = candidate.get("structure_fetch_status", "not_attempted")
            lines.append(
                "| "
                f"{candidate.get('retrieval_rank')} | "
                f"{candidate.get('compound_name')} | "
                f"{candidate.get('design_status')} | "
                f"{candidate.get('evidence_level')} | "
                f"{candidate.get('retrieval_rule_id')} | "
                f"{structure_status} |"
            )
        lines.append("")

    lines.append("## Candidate Reasoning")
    lines.append("")
    for candidate in candidates:
        lines.append(f"### {candidate.get('compound_name')}")
        lines.append("")
        lines.append(f"- Retrieval rank: {candidate.get('retrieval_rank')}")
        lines.append(f"- Retrieval reason: {candidate.get('retrieval_reason')}")
        lines.append(f"- Target family basis: {candidate.get('target_family_basis')}")
        lines.append(f"- Domain basis: {', '.join(candidate.get('domain_basis') or [])}")
        lines.append(f"- Retrieval terms: {', '.join(candidate.get('retrieval_terms') or [])}")
        lines.append(f"- Source databases: {', '.join(candidate.get('source_databases') or [])}")
        if candidate.get("pubchem_cid"):
            lines.append(f"- PubChem CID: {candidate.get('pubchem_cid')}")
        if candidate.get("local_sdf_path"):
            lines.append(f"- Local SDF: {candidate.get('local_sdf_path')}")
        if candidate.get("structure_fetch_error"):
            lines.append(f"- Structure fetch error: {candidate.get('structure_fetch_error')}")
        lines.append("")

    lines.append("## Limitations")
    lines.append("")
    lines.append("- Retrieved compounds are candidates for computational docking, not confirmed antivirals for the submitted target.")
    lines.append("- Stage 4A does not yet evaluate mutation-aware ligand adaptation.")
    lines.append("- Stage 4A does not yet perform pocket-feature-driven novel ligand design.")
    lines.append("- Future versions should add ChEMBL, BindingDB, RCSB ligand extraction, pocket features, and mutation effects.")
    lines.append("")

    return "\n".join(lines)


def run_compound_retrieval(
    *,
    target_evidence_path: Path,
    output_dir: Path,
    homolog_summary_path: Path | None = None,
    max_candidates: int = 20,
    fetch_structures: bool = True,
    pubchem_timeout_seconds: int = 60,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    target_evidence = load_json(target_evidence_path)
    homolog_summary = (
        load_json(homolog_summary_path)
        if homolog_summary_path and homolog_summary_path.exists()
        else None
    )

    candidates = build_ligand_candidates(
        target_evidence,
        homolog_summary,
        max_candidates=max_candidates,
    )

    if fetch_structures and candidates:
        candidates = fetch_candidate_structures_from_pubchem(
            candidates,
            output_dir=output_dir,
            timeout_seconds=pubchem_timeout_seconds,
        )

    candidate_json = output_dir / "ligand_candidates.json"
    candidate_csv = output_dir / "candidate_ligands.csv"
    docking_manifest = output_dir / "docking_manifest.csv"
    report_path = output_dir / "ligand_search_report.md"

    write_json(candidate_json, {"candidates": candidates})
    write_candidate_csv(candidate_csv, candidates)
    write_docking_manifest(docking_manifest, candidates)

    report = render_ligand_search_report(
        target_evidence=target_evidence,
        candidates=candidates,
        fetch_structures=fetch_structures,
    )
    report_path.write_text(report, encoding="utf-8")

    return {
        "ligand_candidates": candidate_json,
        "candidate_ligands_csv": candidate_csv,
        "docking_manifest": docking_manifest,
        "ligand_search_report": report_path,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Stage 4A compound retrieval from target evidence.")
    parser.add_argument("--target-evidence", required=True, type=Path)
    parser.add_argument("--homolog-summary", type=Path, default=None)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-candidates", type=int, default=20)
    parser.add_argument("--pubchem-timeout-seconds", type=int, default=60)
    parser.add_argument(
        "--no-fetch-structures",
        action="store_true",
        help="Generate candidate packets without contacting PubChem for SDF structures.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    outputs = run_compound_retrieval(
        target_evidence_path=args.target_evidence,
        homolog_summary_path=args.homolog_summary,
        output_dir=args.output_dir,
        max_candidates=args.max_candidates,
        fetch_structures=not args.no_fetch_structures,
        pubchem_timeout_seconds=args.pubchem_timeout_seconds,
    )

    print("[COMPOUND_RETRIEVAL] Outputs:")
    for label, path in outputs.items():
        print(f"[COMPOUND_RETRIEVAL] {label}: {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
