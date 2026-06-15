from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any

import numpy as np

from .models import PocketDefinition
from .subprocess_utils import resolve_executable, run_command


def _read_pdb_coordinates(path: Path) -> np.ndarray:
    coordinates: list[tuple[float, float, float]] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            coordinates.append(
                (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            )
        except ValueError:
            continue
    if not coordinates:
        raise RuntimeError(f"No coordinates found in pocket file: {path}")
    return np.asarray(coordinates, dtype=float)


def parse_fpocket_info(info_path: Path) -> list[dict[str, Any]]:
    pockets: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    pocket_pattern = re.compile(r"^Pocket\s+(\d+)\s*:", re.IGNORECASE)
    for raw_line in info_path.read_text(errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        match = pocket_pattern.match(line)
        if match:
            if current:
                pockets.append(current)
            current = {"number": int(match.group(1)), "metrics": {}}
            continue
        if current is not None and ":" in line:
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()
            try:
                parsed: Any = float(value)
            except ValueError:
                parsed = value
            current["metrics"][key] = parsed
    if current:
        pockets.append(current)
    return pockets


def _pocket_score(pocket: dict[str, Any]) -> float:
    metrics = pocket.get("metrics", {})
    for key in ("Score", "Druggability Score", "Drug Score"):
        value = metrics.get(key)
        if isinstance(value, (int, float)):
            return float(value)
    return float("-inf")


def detect_fpocket_box(
    receptor_pdb: Path,
    work_dir: Path,
    *,
    padding: float = 4.0,
    pocket_number: int | None = None,
    fpocket_bin: str = "fpocket",
) -> PocketDefinition:
    work_dir.mkdir(parents=True, exist_ok=True)
    fpocket = resolve_executable(fpocket_bin, "fpocket")
    copied_pdb = work_dir / "fpocket_receptor.pdb"
    shutil.copy2(receptor_pdb, copied_pdb)
    run_command([fpocket, "-f", str(copied_pdb)], cwd=work_dir)

    out_dir = work_dir / f"{copied_pdb.stem}_out"
    info_files = sorted(out_dir.glob("*_info.txt"))
    if not info_files:
        raise RuntimeError("fpocket did not create an *_info.txt file")
    pockets = parse_fpocket_info(info_files[0])
    if not pockets:
        raise RuntimeError("fpocket did not report any pockets")

    if pocket_number is None:
        selected = max(pockets, key=_pocket_score)
        pocket_number = int(selected["number"])
    else:
        matches = [item for item in pockets if item["number"] == pocket_number]
        if not matches:
            raise ValueError(f"fpocket did not report pocket {pocket_number}")
        selected = matches[0]

    candidates = [
        out_dir / "pockets" / f"pocket{pocket_number}_atm.pdb",
        out_dir / "pockets" / f"pocket{pocket_number}_vert.pqr",
        out_dir / "pockets" / f"pocket{pocket_number}_vert.pdb",
    ]
    pocket_file = next((path for path in candidates if path.is_file()), None)
    if pocket_file is None:
        raise RuntimeError(
            f"No coordinate file found for fpocket pocket {pocket_number}"
        )

    coordinates = _read_pdb_coordinates(pocket_file)
    minimum = coordinates.min(axis=0)
    maximum = coordinates.max(axis=0)
    center = (minimum + maximum) / 2.0
    size = maximum - minimum + 2.0 * padding
    size = np.maximum(size, np.asarray([12.0, 12.0, 12.0]))

    return PocketDefinition(
        mode="explicit",
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(size[0]),
        size_y=float(size[1]),
        size_z=float(size[2]),
        source=(
            f"fpocket pocket {pocket_number}; score={_pocket_score(selected):.4f}; "
            f"file={pocket_file}"
        ),
    )


def build_pocket_definition(
    *,
    receptor_pdb: Path,
    work_dir: Path,
    explicit_values: tuple[
        float | None,
        float | None,
        float | None,
        float | None,
        float | None,
        float | None,
    ],
    autobox_ligand: Path | None,
    fpocket_padding: float,
    fpocket_pocket: int | None,
    fpocket_bin: str,
) -> PocketDefinition:
    explicit_count = sum(value is not None for value in explicit_values)
    if autobox_ligand is not None and explicit_count:
        raise ValueError("Use either explicit box values or --autobox-ligand")
    if autobox_ligand is not None:
        return PocketDefinition(
            mode="autobox",
            autobox_ligand=autobox_ligand,
            source=f"autobox ligand {autobox_ligand}",
        )
    if explicit_count:
        if explicit_count != 6:
            raise ValueError(
                "Explicit pocket mode requires center x/y/z and size x/y/z"
            )
        cx, cy, cz, sx, sy, sz = explicit_values
        return PocketDefinition(
            mode="explicit",
            center_x=float(cx),
            center_y=float(cy),
            center_z=float(cz),
            size_x=float(sx),
            size_y=float(sy),
            size_z=float(sz),
            source="user-specified box",
        )
    return detect_fpocket_box(
        receptor_pdb,
        work_dir,
        padding=fpocket_padding,
        pocket_number=fpocket_pocket,
        fpocket_bin=fpocket_bin,
    )
