from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

from folding import run_colabfold
from homolog_search import DEFAULT_API_URL, parse_cpu_response, post_fasta


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

    cpu_dir = output_dir / "cpu"
    colabfold_dir = output_dir / "colabfold"

    print(f"FASTA: {fasta_path}")

    if not args.skip_cpu:
        raw_response = post_fasta(args.api_url, fasta_path)
        parsed_cpu = parse_cpu_response(raw_response)

        raw_path = cpu_dir / f"{fasta_path.stem}_post_response.json"
        parsed_path = cpu_dir / f"{fasta_path.stem}_parsed_cpu_data.json"

        save_json(raw_path, raw_response)
        save_json(parsed_path, parsed_cpu)

        print("\nCPU API complete.")
        print(f"Job ID: {parsed_cpu.get('job_id')}")
        print(f"Status: {parsed_cpu.get('status')}")
        print(f"CDD rows: {parsed_cpu['result_counts']['cdd']}")
        print(f"InterPro rows: {parsed_cpu['result_counts']['interpro']}")
        print(f"VOGDB rows: {parsed_cpu['result_counts']['vogdb']}")
        print(f"Raw CPU response: {raw_path}")
        print(f"Parsed CPU data: {parsed_path}")

    if not args.skip_colabfold:
        best_pdb = run_colabfold(
            fasta_path=fasta_path,
            output_dir=colabfold_dir,
            quick_test=not args.full_colabfold,
            overwrite=args.overwrite_colabfold,
        )

        print("\nColabFold complete.")
        print(f"Best PDB: {best_pdb}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
