"""Utilities for generic protein-ligand reference complex handling.

This module is intentionally target-agnostic. It supports pose-recovery
benchmarks for any known protein-ligand complex, not just HIV protease
or darunavir.

Typical use:
    python -m compoundrank.reference_complex --complex-pdb complex.pdb --list-ligands

    python -m compoundrank.reference_complex \
      --complex-pdb complex.pdb \
      --ligand-resname 017 \
      --output-dir benchmark_extract
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


WATER_RESNAMES = {"HOH", "WAT", "DOD"}


@dataclass(frozen=True)
class LigandKey:
    resname: str
    chain: str
    resseq: str
    icode: str
    altloc: str


@dataclass
class LigandSummary:
    resname: str
    chain: str
    resseq: str
    icode: str
    altloc: str
    atom_count: int
    center_x: float
    center_y: float
    center_z: float
    min_x: float
    min_y: float
    min_z: float
    max_x: float
    max_y: float
    max_z: float


def parse_atom_line(line: str) -> dict[str, Any] | None:
    """Parse ATOM/HETATM line using fixed-width PDB columns."""
    record = line[0:6].strip()
    if record not in {"ATOM", "HETATM"}:
        return None

    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None

    return {
        "record": record,
        "atom_name": line[12:16].strip(),
        "altloc": line[16].strip() or "-",
        "resname": line[17:20].strip(),
        "chain": line[21].strip() or "-",
        "resseq": line[22:26].strip(),
        "icode": line[26].strip() or "-",
        "x": x,
        "y": y,
        "z": z,
        "raw": line.rstrip("\n"),
    }


def ligand_key(atom: dict[str, Any]) -> LigandKey:
    return LigandKey(
        resname=str(atom["resname"]),
        chain=str(atom["chain"]),
        resseq=str(atom["resseq"]),
        icode=str(atom["icode"]),
        altloc=str(atom["altloc"]),
    )


def is_reference_ligand_atom(atom: dict[str, Any]) -> bool:
    if atom["record"] != "HETATM":
        return False
    if str(atom["resname"]).upper() in WATER_RESNAMES:
        return False
    return True


def read_atoms(pdb_path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    for line in pdb_path.read_text(errors="replace").splitlines():
        atom = parse_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    return atoms


def summarize_ligands(pdb_path: Path) -> list[LigandSummary]:
    grouped: dict[LigandKey, list[dict[str, Any]]] = defaultdict(list)

    for atom in read_atoms(pdb_path):
        if is_reference_ligand_atom(atom):
            grouped[ligand_key(atom)].append(atom)

    summaries: list[LigandSummary] = []
    for key, atoms in grouped.items():
        xs = [float(atom["x"]) for atom in atoms]
        ys = [float(atom["y"]) for atom in atoms]
        zs = [float(atom["z"]) for atom in atoms]

        summaries.append(
            LigandSummary(
                resname=key.resname,
                chain=key.chain,
                resseq=key.resseq,
                icode=key.icode,
                altloc=key.altloc,
                atom_count=len(atoms),
                center_x=sum(xs) / len(xs),
                center_y=sum(ys) / len(ys),
                center_z=sum(zs) / len(zs),
                min_x=min(xs),
                min_y=min(ys),
                min_z=min(zs),
                max_x=max(xs),
                max_y=max(ys),
                max_z=max(zs),
            )
        )

    summaries.sort(
        key=lambda item: (
            item.resname,
            item.chain,
            int(item.resseq) if item.resseq.lstrip("-").isdigit() else item.resseq,
            item.altloc,
        )
    )
    return summaries


def select_ligand_key(
    summaries: list[LigandSummary],
    *,
    resname: str,
    chain: str | None = None,
    resseq: str | None = None,
    altloc: str | None = None,
) -> LigandKey:
    matches = []
    for summary in summaries:
        if summary.resname != resname:
            continue
        if chain is not None and summary.chain != chain:
            continue
        if resseq is not None and summary.resseq != resseq:
            continue
        if altloc is not None and summary.altloc != altloc:
            continue
        matches.append(summary)

    if not matches:
        available = "\n".join(
            f"- resname={s.resname} chain={s.chain} resseq={s.resseq} altloc={s.altloc} atoms={s.atom_count}"
            for s in summaries
        )
        raise ValueError(f"No matching ligand found for resname={resname}. Available ligands:\n{available}")

    if len(matches) > 1:
        available = "\n".join(
            f"- resname={s.resname} chain={s.chain} resseq={s.resseq} altloc={s.altloc} atoms={s.atom_count}"
            for s in matches
        )
        raise ValueError(
            "Multiple matching ligands found. Specify --ligand-chain, --ligand-resseq, or --ligand-altloc.\n"
            f"{available}"
        )

    chosen = matches[0]
    return LigandKey(
        resname=chosen.resname,
        chain=chosen.chain,
        resseq=chosen.resseq,
        icode=chosen.icode,
        altloc=chosen.altloc,
    )


def make_box_from_atoms(
    ligand_atoms: list[dict[str, Any]],
    *,
    padding: float,
    min_size: float,
) -> dict[str, Any]:
    xs = [float(atom["x"]) for atom in ligand_atoms]
    ys = [float(atom["y"]) for atom in ligand_atoms]
    zs = [float(atom["z"]) for atom in ligand_atoms]

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    center_x = sum(xs) / len(xs)
    center_y = sum(ys) / len(ys)
    center_z = sum(zs) / len(zs)

    size_x = max(max_x - min_x + 2 * padding, min_size)
    size_y = max(max_y - min_y + 2 * padding, min_size)
    size_z = max(max_z - min_z + 2 * padding, min_size)

    return {
        "box_mode": "reference_ligand",
        "center_x": center_x,
        "center_y": center_y,
        "center_z": center_z,
        "size_x": size_x,
        "size_y": size_y,
        "size_z": size_z,
        "padding": padding,
        "min_size": min_size,
        "ligand_atom_count": len(ligand_atoms),
        "ligand_bounds": {
            "min_x": min_x,
            "min_y": min_y,
            "min_z": min_z,
            "max_x": max_x,
            "max_y": max_y,
            "max_z": max_z,
        },
    }


def extract_reference_complex(
    *,
    complex_pdb: Path,
    output_dir: Path,
    ligand_resname: str,
    ligand_chain: str | None = None,
    ligand_resseq: str | None = None,
    ligand_altloc: str | None = None,
    padding: float = 8.0,
    min_size: float = 18.0,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = read_atoms(complex_pdb)
    summaries = summarize_ligands(complex_pdb)
    selected_key = select_ligand_key(
        summaries,
        resname=ligand_resname,
        chain=ligand_chain,
        resseq=ligand_resseq,
        altloc=ligand_altloc,
    )

    receptor_lines: list[str] = []
    ligand_lines: list[str] = []
    ligand_atoms: list[dict[str, Any]] = []

    for atom in atoms:
        if atom["record"] == "ATOM":
            receptor_lines.append(atom["raw"])
            continue

        if atom["record"] == "HETATM" and is_reference_ligand_atom(atom):
            if ligand_key(atom) == selected_key:
                ligand_lines.append(atom["raw"])
                ligand_atoms.append(atom)

    if not receptor_lines:
        raise ValueError(f"No protein ATOM records found in {complex_pdb}")
    if not ligand_lines:
        raise ValueError(f"No ligand atoms extracted for {selected_key}")

    receptor_path = output_dir / "reference_receptor_protein_only.pdb"
    ligand_path = output_dir / f"reference_ligand_{selected_key.resname}_{selected_key.chain}_{selected_key.resseq}_{selected_key.altloc}.pdb"
    box_path = output_dir / "reference_box.json"
    summary_path = output_dir / "reference_ligand_summary.json"

    receptor_path.write_text("\n".join(receptor_lines) + "\nEND\n", encoding="utf-8")
    ligand_path.write_text("\n".join(ligand_lines) + "\nEND\n", encoding="utf-8")

    box = make_box_from_atoms(ligand_atoms, padding=padding, min_size=min_size)
    box["reference_complex"] = str(complex_pdb)
    box["reference_ligand"] = asdict(selected_key)
    box_path.write_text(json.dumps(box, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    summary = [asdict(item) for item in summaries]
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "reference_receptor": receptor_path,
        "reference_ligand": ligand_path,
        "reference_box": box_path,
        "reference_ligand_summary": summary_path,
    }


def print_ligand_table(summaries: list[LigandSummary]) -> None:
    print("resname,chain,resseq,icode,altloc,atom_count,center_x,center_y,center_z")
    for item in summaries:
        print(
            f"{item.resname},{item.chain},{item.resseq},{item.icode},{item.altloc},"
            f"{item.atom_count},{item.center_x:.3f},{item.center_y:.3f},{item.center_z:.3f}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="List/extract ligands from a reference protein-ligand PDB complex.")
    parser.add_argument("--complex-pdb", required=True, type=Path)
    parser.add_argument("--list-ligands", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--ligand-resname")
    parser.add_argument("--ligand-chain")
    parser.add_argument("--ligand-resseq")
    parser.add_argument("--ligand-altloc")
    parser.add_argument("--padding", type=float, default=8.0)
    parser.add_argument("--min-size", type=float, default=18.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    summaries = summarize_ligands(args.complex_pdb)

    if args.list_ligands:
        print_ligand_table(summaries)
        return 0

    if not args.output_dir:
        raise SystemExit("--output-dir is required unless --list-ligands is used")
    if not args.ligand_resname:
        raise SystemExit("--ligand-resname is required unless --list-ligands is used")

    outputs = extract_reference_complex(
        complex_pdb=args.complex_pdb,
        output_dir=args.output_dir,
        ligand_resname=args.ligand_resname,
        ligand_chain=args.ligand_chain,
        ligand_resseq=args.ligand_resseq,
        ligand_altloc=args.ligand_altloc,
        padding=args.padding,
        min_size=args.min_size,
    )

    print("[REFERENCE_COMPLEX] Outputs:")
    for label, path in outputs.items():
        print(f"[REFERENCE_COMPLEX] {label}: {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
