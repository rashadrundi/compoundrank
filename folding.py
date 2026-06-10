from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def run_colabfold(
    fasta_path: Path,
    output_dir: Path,
    quick_test: bool = True,
    overwrite: bool = False,
) -> Path:
    if shutil.which("colabfold_batch") is None:
        raise RuntimeError("colabfold_batch was not found on PATH")

    output_dir.mkdir(parents=True, exist_ok=True)

    existing_pdbs = sorted(output_dir.glob("*.pdb"))
    if existing_pdbs and not overwrite:
        rank_001 = [p for p in existing_pdbs if "rank_001" in p.name]
        best_existing = rank_001[0] if rank_001 else existing_pdbs[0]
        print(f"Using existing ColabFold PDB: {best_existing}")
        return best_existing

    command = ["colabfold_batch"]

    if quick_test:
        command += [
            "--msa-mode", "single_sequence",
            "--num-models", "1",
            "--num-recycle", "1",
        ]

    if overwrite:
        command.append("--overwrite-existing-results")

    command += [
        str(fasta_path),
        str(output_dir),
    ]

    print("Running ColabFold:")
    print(" ".join(command))

    subprocess.run(command, check=True)

    pdb_files = sorted(output_dir.glob("*.pdb"))

    if not pdb_files:
        raise RuntimeError(f"ColabFold finished, but no PDB files were found in {output_dir}")

    rank_001 = [p for p in pdb_files if "rank_001" in p.name]
    return rank_001[0] if rank_001 else pdb_files[0]
