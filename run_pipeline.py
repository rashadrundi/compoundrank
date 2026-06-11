from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from folding import run_colabfold
from homolog_search import DEFAULT_API_URL, parse_cpu_response, post_fasta
from structure_quality import write_quality_report
from pocket_detection import run_fpocket


FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna"}


def find_single_fasta(input_dir: Path) -> Path:
    fasta_files = [
        path for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    ]

    if not fasta_files:
        raise FileNotFoundError(f"No FASTA file found in {input_dir}")

    if len(fasta_files) > 1:
        names = "\n".join(f"  - {p}" for p in fasta_files)
        raise ValueError(f"More than one FASTA found. Use --fasta explicitly.\n{names}")

    return fasta_files[0]


def save_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")

def run_cpu_branch(
    api_url: str,
    fasta_path: Path,
    cpu_dir: Path,
) -> dict[str, Any]:
    print("\n[CPU] Sending FASTA to CPU API...")

    raw_response = post_fasta(api_url, fasta_path)
    parsed_cpu = parse_cpu_response(raw_response)

    raw_path = cpu_dir / f"{fasta_path.stem}_post_response.json"
    parsed_path = cpu_dir / f"{fasta_path.stem}_parsed_cpu_data.json"

    save_json(raw_path, raw_response)
    save_json(parsed_path, parsed_cpu)

    print("\n[CPU] API complete.")
    print(f"[CPU] Job ID: {parsed_cpu.get('job_id')}")
    print(f"[CPU] Status: {parsed_cpu.get('status')}")
    print(f"[CPU] CDD rows: {parsed_cpu['result_counts']['cdd']}")
    print(f"[CPU] InterPro rows: {parsed_cpu['result_counts']['interpro']}")
    print(f"[CPU] VOGDB rows: {parsed_cpu['result_counts']['vogdb']}")
    print(f"[CPU] Raw response: {raw_path}")
    print(f"[CPU] Parsed data: {parsed_path}")

    return parsed_cpu


def run_structure_branch(
    fasta_path: Path,
    output_dir: Path,
    quick_test: bool,
    overwrite_colabfold: bool,
) -> dict[str, Any]:
    colabfold_dir = output_dir / "colabfold"
    quality_path = output_dir / "quality" / "structure_quality.json"
    pocket_path = output_dir / "pockets" / "fpocket_summary.json"

    print("\n[STRUCTURE] Starting local structure branch...")

    best_pdb = run_colabfold(
        fasta_path=fasta_path,
        output_dir=colabfold_dir,
        quick_test=quick_test,
        overwrite=overwrite_colabfold,
    )

    print("\n[STRUCTURE] ColabFold complete.")
    print(f"[STRUCTURE] Best PDB: {best_pdb}")

    quality_report = write_quality_report(
        colabfold_dir=colabfold_dir,
        output_path=quality_path,
        pdb_path=best_pdb,
    )

    print("\n[STRUCTURE] Quality extraction complete.")
    print(
        f"[STRUCTURE] Mean pLDDT: "
        f"{quality_report['plddt']['mean_plddt']}"
    )
    print(
        f"[STRUCTURE] PAE available: "
        f"{quality_report['pae']['available']}"
    )
    print(f"[STRUCTURE] Quality report: {quality_path}")

    pocket_result = run_fpocket(
        pdb_path=best_pdb,
        output_json=pocket_path,
        overwrite=False,
    )

    print("\n[STRUCTURE] Pocket detection complete.")
    print(f"[STRUCTURE] Pocket count: {pocket_result['pocket_count']}")
    print(f"[STRUCTURE] Pocket summary: {pocket_path}")

    return {
        "best_pdb": str(best_pdb),
        "quality": quality_report,
        "pockets": pocket_result,
    }

def main() -> int:
    parser = argparse.ArgumentParser(
        description="CompoundRank local runner: CPU homolog search + local ColabFold."
    )

    parser.add_argument(
        "--api-url",
        default=os.environ.get("COMPOUNDRANK_API_URL", DEFAULT_API_URL),
    )

    parser.add_argument(
        "--input-dir",
        default="input",
    )

    parser.add_argument(
        "--fasta",
        default=None,
    )

    parser.add_argument(
        "--output-dir",
        default="output",
    )

    parser.add_argument(
        "--skip-cpu",
        action="store_true",
        help="Skip CPU API POST request.",
    )

    parser.add_argument(
        "--skip-colabfold",
        action="store_true",
        help="Skip local ColabFold prediction.",
    )

    parser.add_argument(
        "--full-colabfold",
        action="store_true",
        help="Use normal ColabFold settings instead of quick smoke-test settings.",
    )

    parser.add_argument(
        "--overwrite-colabfold",
        action="store_true",
        help="Force ColabFold to recompute even if a PDB already exists.",
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    fasta_path = Path(args.fasta) if args.fasta else find_single_fasta(input_dir)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file does not exist: {fasta_path}")

    print(f"FASTA: {fasta_path}")

    cpu_dir = output_dir / "cpu"

    futures = {}
    results: dict[str, Any] = {}
    errors: dict[str, Exception] = {}

    print("\nStarting available pipeline branches concurrently...")

    with ThreadPoolExecutor(max_workers=2) as executor:
        if not args.skip_cpu:
            cpu_future = executor.submit(
                run_cpu_branch,
                args.api_url,
                fasta_path,
                cpu_dir,
            )
            futures[cpu_future] = "cpu"

        if not args.skip_colabfold:
            structure_future = executor.submit(
                run_structure_branch,
                fasta_path,
                output_dir,
                not args.full_colabfold,
                args.overwrite_colabfold,
            )
            futures[structure_future] = "structure"

        for future in as_completed(futures):
            branch_name = futures[future]

            try:
                results[branch_name] = future.result()
                print(f"\n[{branch_name.upper()}] Branch finished successfully.")

            except Exception as error:
                errors[branch_name] = error
                print(
                    f"\n[{branch_name.upper()}] Branch failed: {error}",
                    file=sys.stderr,
                )

    if errors:
        details = "; ".join(
            f"{name}: {error}"
            for name, error in errors.items()
        )
        raise RuntimeError(f"One or more pipeline branches failed: {details}")

    print("\nPipeline complete.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
