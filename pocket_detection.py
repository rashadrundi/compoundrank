from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any


def find_best_pdb(colabfold_dir: Path) -> Path:
    pdbs = sorted(colabfold_dir.rglob("*.pdb"))

    # Avoid fpocket-generated pocket PDBs.
    pdbs = [
        p for p in pdbs
        if "_out" not in str(p)
        and "pocket" not in p.name.lower()
    ]

    if not pdbs:
        raise FileNotFoundError(f"No source PDB files found in {colabfold_dir}")

    rank_001 = [p for p in pdbs if "rank_001" in p.name]
    return rank_001[0] if rank_001 else pdbs[0]


def expected_fpocket_out_dir(pdb_path: Path) -> Path:
    return pdb_path.with_name(f"{pdb_path.stem}_out")


def parse_fpocket_info(info_path: Path) -> list[dict[str, Any]]:
    pockets: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    for raw_line in info_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if line.lower().startswith("pocket") and ":" in line:
            if current:
                pockets.append(current)

            pocket_name = line.split(":", 1)[0].strip()
            current = {
                "name": pocket_name,
                "metrics": {},
            }
            continue

        if current is not None and ":" in line:
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()

            try:
                parsed_value: Any = float(value)
                if parsed_value.is_integer():
                    parsed_value = int(parsed_value)
            except ValueError:
                parsed_value = value

            current["metrics"][key] = parsed_value

    if current:
        pockets.append(current)

    return pockets


def run_fpocket(
    pdb_path: Path,
    output_json: Path,
    overwrite: bool = False,
) -> dict[str, Any]:
    if shutil.which("fpocket") is None:
        raise RuntimeError("fpocket was not found on PATH")

    pdb_path = pdb_path.resolve()
    out_dir = expected_fpocket_out_dir(pdb_path)

    if out_dir.exists() and not overwrite:
        print(f"Using existing fpocket output: {out_dir}")
    else:
        if out_dir.exists() and overwrite:
            shutil.rmtree(out_dir)

        command = ["fpocket", "-f", str(pdb_path)]

        print("Running fpocket:")
        print(" ".join(command))

        subprocess.run(command, check=True)

    info_files = sorted(out_dir.glob("*_info.txt"))
    if not info_files:
        raise RuntimeError(f"fpocket finished, but no *_info.txt file was found in {out_dir}")

    info_path = info_files[0]
    pockets = parse_fpocket_info(info_path)

    pocket_files = sorted((out_dir / "pockets").glob("*")) if (out_dir / "pockets").exists() else []

    result = {
        "source_pdb": str(pdb_path),
        "fpocket_output_dir": str(out_dir),
        "info_file": str(info_path),
        "pocket_count": len(pockets),
        "pockets": pockets,
        "pocket_files": [str(p) for p in pocket_files],
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2), encoding="utf-8")

    return result
