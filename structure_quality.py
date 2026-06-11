from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import mean
from typing import Any


def find_best_pdb(colabfold_dir: Path) -> Path:
    pdbs = sorted(colabfold_dir.rglob("*.pdb"))

    if not pdbs:
        raise FileNotFoundError(f"No PDB files found in {colabfold_dir}")

    rank_001 = [p for p in pdbs if "rank_001" in p.name]
    return rank_001[0] if rank_001 else pdbs[0]


def parse_plddt_from_pdb(pdb_path: Path) -> dict[str, Any]:
    residue_scores: dict[str, float] = {}

    with pdb_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            residue_name = line[17:20].strip()
            chain_id = line[21].strip() or "_"
            residue_number = line[22:26].strip()

            try:
                plddt = float(line[60:66].strip())
            except ValueError:
                continue

            key = f"{chain_id}:{residue_number}:{residue_name}"
            residue_scores[key] = plddt

    if not residue_scores:
        raise RuntimeError(f"No CA atom pLDDT values found in {pdb_path}")

    scores = list(residue_scores.values())

    return {
        "source": "pdb_b_factor_column",
        "pdb_file": str(pdb_path),
        "residue_count": len(scores),
        "mean_plddt": round(mean(scores), 2),
        "min_plddt": round(min(scores), 2),
        "max_plddt": round(max(scores), 2),
        "counts": {
            "very_high_90_plus": sum(v >= 90 for v in scores),
            "confident_70_to_90": sum(70 <= v < 90 for v in scores),
            "low_50_to_70": sum(50 <= v < 70 for v in scores),
            "very_low_below_50": sum(v < 50 for v in scores),
        },
        "low_confidence_residues_below_70": {
            k: v for k, v in residue_scores.items() if v < 70
        },
    }


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_matrix(value: Any) -> bool:
    if not isinstance(value, list) or not value:
        return False

    first_row = value[0]
    if not isinstance(first_row, list) or not first_row:
        return False

    first_value = first_row[0]
    return isinstance(first_value, (int, float))


def find_first_key_recursive(obj: Any, key_names: set[str]) -> Any | None:
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key in key_names:
                return value

        for value in obj.values():
            found = find_first_key_recursive(value, key_names)
            if found is not None:
                return found

    elif isinstance(obj, list):
        for value in obj:
            found = find_first_key_recursive(value, key_names)
            if found is not None:
                return found

    return None


def find_score_json(colabfold_dir: Path) -> Path | None:
    json_files = sorted(colabfold_dir.rglob("*.json"))

    if not json_files:
        return None

    preferred = [
        p for p in json_files
        if "score" in p.name.lower() or "pae" in p.name.lower()
    ]

    candidates = preferred or json_files

    for path in candidates:
        try:
            data = load_json(path)
        except Exception:
            continue

        pae_value = find_first_key_recursive(
            data,
            {"pae", "predicted_aligned_error", "predicted_aligned_errors"},
        )

        plddt_value = find_first_key_recursive(data, {"plddt", "max_pae", "ptm"})

        if pae_value is not None or plddt_value is not None:
            return path

    return candidates[0]


def summarize_pae(score_json_path: Path | None) -> dict[str, Any]:
    if score_json_path is None:
        return {
            "available": False,
            "reason": "No ColabFold JSON files found.",
        }

    try:
        data = load_json(score_json_path)
    except Exception as error:
        return {
            "available": False,
            "score_json_file": str(score_json_path),
            "reason": f"Could not parse JSON: {error}",
        }

    pae_value = find_first_key_recursive(
        data,
        {"pae", "predicted_aligned_error", "predicted_aligned_errors"},
    )

    ptm = find_first_key_recursive(data, {"ptm"})
    max_pae = find_first_key_recursive(data, {"max_pae"})

    result: dict[str, Any] = {
        "available": False,
        "score_json_file": str(score_json_path),
        "ptm": ptm,
        "max_pae": max_pae,
    }

    if not is_numeric_matrix(pae_value):
        result["reason"] = "No PAE matrix found in this JSON."
        return result

    flat = [float(x) for row in pae_value for x in row]

    result.update({
        "available": True,
        "matrix_size": {
            "rows": len(pae_value),
            "cols": len(pae_value[0]) if pae_value else 0,
        },
        "mean_pae": round(mean(flat), 2),
        "min_pae": round(min(flat), 2),
        "max_pae_observed": round(max(flat), 2),
        "fractions": {
            "pae_below_5": round(sum(v < 5 for v in flat) / len(flat), 3),
            "pae_below_10": round(sum(v < 10 for v in flat) / len(flat), 3),
            "pae_below_20": round(sum(v < 20 for v in flat) / len(flat), 3),
        },
    })

    return result


def write_quality_report(
    colabfold_dir: Path,
    output_path: Path,
    pdb_path: Path | None = None,
) -> dict[str, Any]:
    pdb = pdb_path if pdb_path is not None else find_best_pdb(colabfold_dir)
    score_json = find_score_json(colabfold_dir)

    report = {
        "pdb_file": str(pdb),
        "colabfold_dir": str(colabfold_dir),
        "plddt": parse_plddt_from_pdb(pdb),
        "pae": summarize_pae(score_json),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract pLDDT and PAE quality information from existing ColabFold output."
    )

    parser.add_argument(
        "--colabfold-dir",
        default="output/colabfold",
    )

    parser.add_argument(
        "--pdb",
        default=None,
    )

    parser.add_argument(
        "--output",
        default="output/quality/structure_quality.json",
    )

    args = parser.parse_args()

    colabfold_dir = Path(args.colabfold_dir)
    pdb_path = Path(args.pdb) if args.pdb else None
    output_path = Path(args.output)

    report = write_quality_report(
        colabfold_dir=colabfold_dir,
        output_path=output_path,
        pdb_path=pdb_path,
    )

    print("Structure confidence extraction complete.")
    print(f"PDB: {report['pdb_file']}")
    print(f"Mean pLDDT: {report['plddt']['mean_plddt']}")
    print(f"Residues checked: {report['plddt']['residue_count']}")
    print(f"PAE available: {report['pae']['available']}")

    if report["pae"]["available"]:
        print(f"Mean PAE: {report['pae']['mean_pae']}")
        print(f"PAE matrix: {report['pae']['matrix_size']}")

    print(f"Saved: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
