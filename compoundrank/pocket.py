from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .models import PocketDefinition
from .paths import sanitize_name
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


def _coordinate_file_for_pocket(out_dir: Path, pocket_number: int) -> Path:
    # The fpocket alpha-sphere vertices describe the pocket
    # geometry directly. Prefer them over nearby receptor atoms so
    # box construction matches the validated blind-redocking workflow.
    candidates = [
        out_dir / "pockets" / f"pocket{pocket_number}_vert.pqr",
        out_dir / "pockets" / f"pocket{pocket_number}_vert.pdb",
        out_dir / "pockets" / f"pocket{pocket_number}_atm.pdb",
    ]
    pocket_file = next((path for path in candidates if path.is_file()), None)
    if pocket_file is None:
        raise RuntimeError(
            f"No coordinate file found for fpocket pocket {pocket_number}"
        )
    return pocket_file


def _box_from_pocket_file(
    pocket_file: Path,
    *,
    padding: float,
    minimum_size: float = 20.0,
) -> tuple[np.ndarray, np.ndarray]:
    if padding < 0:
        raise ValueError("Pocket padding cannot be negative")

    if minimum_size <= 0:
        raise ValueError(
            "Pocket minimum size must be greater than zero"
        )

    coordinates = _read_pdb_coordinates(pocket_file)
    minimum = coordinates.min(axis=0)
    maximum = coordinates.max(axis=0)
    center = (minimum + maximum) / 2.0
    size = maximum - minimum + 2.0 * padding

    size = np.maximum(
        size,
        np.asarray(
            [minimum_size, minimum_size, minimum_size],
            dtype=float,
        ),
    )

    return center, size



def pocket_definition_to_dict(
    pocket: PocketDefinition,
) -> dict[str, Any]:
    return {
        "pocket_id": pocket.pocket_id,
        "pocket_rank": pocket.pocket_rank,
        "mode": pocket.mode,
        "source": pocket.source,
        "fpocket_score": pocket.fpocket_score,
        "center_x": pocket.center_x,
        "center_y": pocket.center_y,
        "center_z": pocket.center_z,
        "size_x": pocket.size_x,
        "size_y": pocket.size_y,
        "size_z": pocket.size_z,
        "autobox_ligand": (
            str(pocket.autobox_ligand)
            if pocket.autobox_ligand is not None
            else None
        ),
    }


def write_pocket_definitions(
    output_path: Path,
    pockets: Iterable[PocketDefinition],
) -> Path:
    pocket_list = list(pockets)

    payload = {
        "pocket_count": len(pocket_list),
        "ranking_method": (
            "fpocket score descending when fpocket is used"
        ),
        "reference_ligand_used_for_selection": False,
        "box_geometry": {
            "coordinate_source_preference": [
                "fpocket vertex PQR",
                "fpocket vertex PDB",
                "fpocket atom PDB",
            ],
            "default_padding_angstrom_per_side": 4.0,
            "minimum_dimension_angstrom": 20.0,
        },
        "pockets": [
            pocket_definition_to_dict(pocket)
            for pocket in pocket_list
        ],
    }

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return output_path


def detect_fpocket_boxes(
    receptor_pdb: Path,
    work_dir: Path,
    *,
    padding: float = 4.0,
    pocket_number: int | None = None,
    top_n: int = 1,
    fpocket_bin: str = "fpocket",
) -> list[PocketDefinition]:
    if top_n < 1:
        raise ValueError("fpocket top_n must be at least 1")

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

    if pocket_number is not None:
        selected_pockets = [
            item for item in pockets if int(item["number"]) == pocket_number
        ]
        if not selected_pockets:
            raise ValueError(f"fpocket did not report pocket {pocket_number}")
    else:
        selected_pockets = sorted(pockets, key=_pocket_score, reverse=True)[:top_n]

    definitions: list[PocketDefinition] = []
    for rank, selected in enumerate(selected_pockets, start=1):
        selected_number = int(selected["number"])
        selected_score = _pocket_score(selected)
        pocket_file = _coordinate_file_for_pocket(out_dir, selected_number)
        center, size = _box_from_pocket_file(pocket_file, padding=padding)
        pocket_id = sanitize_name(f"fpocket_{rank:02d}_pocket_{selected_number}")

        definitions.append(
            PocketDefinition(
                mode="explicit",
                center_x=float(center[0]),
                center_y=float(center[1]),
                center_z=float(center[2]),
                size_x=float(size[0]),
                size_y=float(size[1]),
                size_z=float(size[2]),
                source=(
                    f"fpocket rank {rank}; pocket {selected_number}; "
                    f"score={selected_score:.4f}"
                ),
                pocket_id=pocket_id,
                pocket_rank=rank,
                fpocket_score=selected_score,
            )
        )

    return definitions


def detect_fpocket_box(
    receptor_pdb: Path,
    work_dir: Path,
    *,
    padding: float = 4.0,
    pocket_number: int | None = None,
    fpocket_bin: str = "fpocket",
) -> PocketDefinition:
    return detect_fpocket_boxes(
        receptor_pdb,
        work_dir,
        padding=padding,
        pocket_number=pocket_number,
        top_n=1,
        fpocket_bin=fpocket_bin,
    )[0]


def build_pocket_definitions(
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
    fpocket_top_n: int,
    fpocket_bin: str,
) -> list[PocketDefinition]:
    explicit_count = sum(value is not None for value in explicit_values)

    if autobox_ligand is not None and explicit_count:
        raise ValueError("Use either explicit box values or --autobox-ligand")

    if autobox_ligand is not None:
        return [
            PocketDefinition(
                mode="autobox",
                autobox_ligand=autobox_ligand,
                source=f"autobox ligand {autobox_ligand}",
                pocket_id="autobox_01",
                pocket_rank=1,
            )
        ]

    if explicit_count:
        if explicit_count != 6:
            raise ValueError(
                "Explicit pocket mode requires center x/y/z and size x/y/z"
            )
        cx, cy, cz, sx, sy, sz = explicit_values
        return [
            PocketDefinition(
                mode="explicit",
                center_x=float(cx),
                center_y=float(cy),
                center_z=float(cz),
                size_x=float(sx),
                size_y=float(sy),
                size_z=float(sz),
                source="user-specified box",
                pocket_id="user_box_01",
                pocket_rank=1,
            )
        ]

    return detect_fpocket_boxes(
        receptor_pdb,
        work_dir,
        padding=fpocket_padding,
        pocket_number=fpocket_pocket,
        top_n=fpocket_top_n,
        fpocket_bin=fpocket_bin,
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
    return build_pocket_definitions(
        receptor_pdb=receptor_pdb,
        work_dir=work_dir,
        explicit_values=explicit_values,
        autobox_ligand=autobox_ligand,
        fpocket_padding=fpocket_padding,
        fpocket_pocket=fpocket_pocket,
        fpocket_top_n=1,
        fpocket_bin=fpocket_bin,
    )[0]
