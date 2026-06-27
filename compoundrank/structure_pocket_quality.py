from __future__ import annotations

import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import gemmi


SCHEMA_VERSION = "structure_pocket_quality.v0.1"
DEFAULT_NEAR_BOX_THRESHOLD_ANGSTROM = 4.0


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)

    data = json.loads(
        path.read_text(encoding="utf-8")
    )

    if not isinstance(data, dict):
        raise TypeError(
            f"Expected JSON object in {path}"
        )

    return data


def _parse_residue_number(
    value: object,
) -> tuple[int, str]:
    match = re.fullmatch(
        r"\s*(-?\d+)([A-Za-z]?)\s*",
        str(value),
    )

    if not match:
        raise ValueError(
            "Unsupported residue number: "
            f"{value!r}"
        )

    return int(match.group(1)), match.group(2)


def _residue_atom_coordinates(
    structure: gemmi.Structure,
    *,
    model_index: int,
    chain_id: str,
    residue_number: object,
) -> list[tuple[float, float, float]]:
    if model_index < 0 or model_index >= len(structure):
        raise IndexError(
            "Model index is outside the "
            f"structure: {model_index}"
        )

    sequence_number, insertion_code = (
        _parse_residue_number(
            residue_number
        )
    )

    coordinates: list[
        tuple[float, float, float]
    ] = []

    model = structure[model_index]

    for chain in model:
        if chain.name != chain_id:
            continue

        for residue in chain:
            if residue.seqid.num != sequence_number:
                continue

            current_insertion_code = (
                residue.seqid.icode.strip()
            )

            if (
                current_insertion_code
                != insertion_code
            ):
                continue

            for atom in residue:
                if atom.element.name == "H":
                    continue

                coordinates.append(
                    (
                        float(atom.pos.x),
                        float(atom.pos.y),
                        float(atom.pos.z),
                    )
                )

    return coordinates


def _explicit_box_geometry(
    pocket: dict[str, Any],
) -> tuple[
    tuple[float, float, float],
    tuple[float, float, float],
]:
    if pocket.get("mode") != "explicit":
        raise ValueError(
            "Pocket-local structural analysis "
            "currently requires an explicit "
            "docking box."
        )

    coordinate_keys = (
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
    )

    missing = [
        key
        for key in coordinate_keys
        if pocket.get(key) is None
    ]

    if missing:
        raise ValueError(
            "Pocket definition is missing "
            f"geometry fields: {missing}"
        )

    center = (
        float(pocket["center_x"]),
        float(pocket["center_y"]),
        float(pocket["center_z"]),
    )

    half_sizes = (
        float(pocket["size_x"]) / 2.0,
        float(pocket["size_y"]) / 2.0,
        float(pocket["size_z"]) / 2.0,
    )

    return center, half_sizes


def _point_to_box_distance(
    point: tuple[float, float, float],
    pocket: dict[str, Any],
) -> float:
    center, half_sizes = (
        _explicit_box_geometry(pocket)
    )

    outside_components = [
        max(
            abs(coordinate - box_center)
            - half_size,
            0.0,
        )
        for coordinate, box_center, half_size
        in zip(
            point,
            center,
            half_sizes,
        )
    ]

    return math.sqrt(
        sum(
            component**2
            for component in outside_components
        )
    )


def _point_to_box_center_distance(
    point: tuple[float, float, float],
    pocket: dict[str, Any],
) -> float:
    center, _ = _explicit_box_geometry(
        pocket
    )

    return math.sqrt(
        sum(
            (
                coordinate
                - box_center
            )
            ** 2
            for coordinate, box_center
            in zip(
                point,
                center,
            )
        )
    )


def _classify_box_distance(
    distance_angstrom: float,
    *,
    near_threshold_angstrom: float,
) -> str:
    if distance_angstrom <= 1e-9:
        return "inside_docking_box"

    if (
        distance_angstrom
        <= near_threshold_angstrom
    ):
        return "near_docking_box"

    return "distal_from_docking_box"


def evaluate_structure_pocket_quality(
    *,
    structure_path: Path,
    ramachandran_report_path: Path,
    pocket_definitions_path: Path,
    pocket_selection_summary_path: Path,
    near_threshold_angstrom: float = (
        DEFAULT_NEAR_BOX_THRESHOLD_ANGSTROM
    ),
) -> dict[str, Any]:
    if near_threshold_angstrom < 0:
        raise ValueError(
            "near_threshold_angstrom cannot "
            "be negative."
        )

    ramachandran = _read_json(
        ramachandran_report_path
    )
    pocket_data = _read_json(
        pocket_definitions_path
    )
    selection_data = _read_json(
        pocket_selection_summary_path
    )

    structure = gemmi.read_structure(
        str(structure_path)
    )

    model_index = int(
        ramachandran.get(
            "model_index",
            0,
        )
        or 0
    )

    outliers = [
        residue
        for residue
        in ramachandran.get(
            "residues",
            [],
        )
        if (
            isinstance(residue, dict)
            and residue.get(
                "classification"
            )
            == "outlier"
        )
    ]

    pockets = {
        str(pocket["pocket_id"]): pocket
        for pocket
        in pocket_data.get(
            "pockets",
            [],
        )
        if (
            isinstance(pocket, dict)
            and pocket.get("pocket_id")
        )
    }

    selected_by_pocket: dict[
        str,
        list[str],
    ] = defaultdict(list)

    for selection in selection_data.get(
        "selected_pockets",
        [],
    ):
        if not isinstance(
            selection,
            dict,
        ):
            continue

        if not selection.get("selected"):
            continue

        pocket_id = str(
            selection["pocket_id"]
        )

        selected_by_pocket[
            pocket_id
        ].append(
            str(selection["compound"])
        )

    selected_pocket_ids = sorted(
        selected_by_pocket
    )

    missing_selected_pocket_ids = [
        pocket_id
        for pocket_id
        in selected_pocket_ids
        if pocket_id not in pockets
    ]

    outlier_results: list[
        dict[str, Any]
    ] = []

    unresolved_outliers: list[str] = []
    inside_selected_box: list[str] = []
    near_selected_box: list[str] = []

    for outlier in outliers:
        residue_label = str(
            outlier["residue"]
        )

        coordinates = (
            _residue_atom_coordinates(
                structure,
                model_index=model_index,
                chain_id=str(
                    outlier["chain"]
                ),
                residue_number=outlier[
                    "residue_number"
                ],
            )
        )

        if not coordinates:
            unresolved_outliers.append(
                residue_label
            )

        pocket_results: list[
            dict[str, Any]
        ] = []

        for pocket_id, pocket in (
            pockets.items()
        ):
            selected_compounds = sorted(
                selected_by_pocket.get(
                    pocket_id,
                    [],
                )
            )

            try:
                if coordinates:
                    minimum_box_distance = min(
                        _point_to_box_distance(
                            coordinate,
                            pocket,
                        )
                        for coordinate
                        in coordinates
                    )

                    minimum_center_distance = (
                        min(
                            _point_to_box_center_distance(
                                coordinate,
                                pocket,
                            )
                            for coordinate
                            in coordinates
                        )
                    )

                    localization = (
                        _classify_box_distance(
                            minimum_box_distance,
                            near_threshold_angstrom=(
                                near_threshold_angstrom
                            ),
                        )
                    )
                else:
                    minimum_box_distance = None
                    minimum_center_distance = None
                    localization = (
                        "residue_not_found"
                    )

                geometry_status = "complete"
                geometry_error = None

            except ValueError as exc:
                minimum_box_distance = None
                minimum_center_distance = None
                localization = (
                    "box_geometry_unavailable"
                )
                geometry_status = "unavailable"
                geometry_error = str(exc)

            pocket_results.append(
                {
                    "pocket_id": pocket_id,
                    "pocket_rank": pocket.get(
                        "pocket_rank"
                    ),
                    "fpocket_score": (
                        pocket.get(
                            "fpocket_score"
                        )
                    ),
                    "selected": bool(
                        selected_compounds
                    ),
                    "selected_compounds": (
                        selected_compounds
                    ),
                    "minimum_distance_to_box_angstrom": (
                        minimum_box_distance
                    ),
                    "minimum_distance_to_center_angstrom": (
                        minimum_center_distance
                    ),
                    "localization": localization,
                    "geometry_status": (
                        geometry_status
                    ),
                    "geometry_error": (
                        geometry_error
                    ),
                }
            )

        selected_results = [
            result
            for result in pocket_results
            if result["selected"]
        ]

        measurable_selected_results = [
            result
            for result in selected_results
            if result[
                "minimum_distance_to_box_angstrom"
            ]
            is not None
        ]

        if measurable_selected_results:
            nearest_selected = min(
                measurable_selected_results,
                key=lambda result: result[
                    "minimum_distance_to_box_angstrom"
                ],
            )

            nearest_localization = (
                nearest_selected[
                    "localization"
                ]
            )

            if (
                nearest_localization
                == "inside_docking_box"
            ):
                inside_selected_box.append(
                    residue_label
                )

            elif (
                nearest_localization
                == "near_docking_box"
            ):
                near_selected_box.append(
                    residue_label
                )
        else:
            nearest_selected = None

        outlier_results.append(
            {
                "residue": residue_label,
                "chain": outlier["chain"],
                "residue_name": outlier[
                    "residue_name"
                ],
                "residue_number": outlier[
                    "residue_number"
                ],
                "ramachandran_category": (
                    outlier["category"]
                ),
                "ramachandran_score": (
                    outlier["score"]
                ),
                "atom_count": len(
                    coordinates
                ),
                "nearest_selected_pocket": (
                    nearest_selected
                ),
                "pockets": pocket_results,
            }
        )

    global_summary = ramachandran.get(
        "summary",
        {},
    )

    if not isinstance(
        global_summary,
        dict,
    ):
        global_summary = {}

    global_goals_met = bool(
        global_summary.get(
            "favored_goal_met"
        )
        and global_summary.get(
            "outlier_goal_met"
        )
    )

    selected_box_local_outliers = sorted(
        set(
            inside_selected_box
            + near_selected_box
        )
    )

    if missing_selected_pocket_ids:
        verdict = (
            "manual_review_missing_"
            "selected_pocket_geometry"
        )

    elif unresolved_outliers:
        verdict = (
            "manual_review_unresolved_"
            "outlier_residues"
        )

    elif inside_selected_box:
        verdict = (
            "selected_pocket_geometry_"
            "concern"
        )

    elif near_selected_box:
        verdict = (
            "manual_review_of_"
            "selected_pocket"
        )

    elif global_goals_met:
        verdict = "strong"

    else:
        verdict = (
            "usable_with_global_"
            "geometry_caution"
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "source_structure": str(
            structure_path
        ),
        "ramachandran_report": str(
            ramachandran_report_path
        ),
        "pocket_definitions": str(
            pocket_definitions_path
        ),
        "pocket_selection_summary": str(
            pocket_selection_summary_path
        ),
        "model_index": model_index,
        "near_box_threshold_angstrom": (
            near_threshold_angstrom
        ),
        "interpretation_scope": (
            "Distances are measured to padded "
            "GNINA docking boxes, not directly "
            "to the molecular pocket surface."
        ),
        "limitations": [
            (
                "A docking box includes empty "
                "search volume and is broader "
                "than the physical pocket."
            ),
            (
                "Distal Ramachandran outliers "
                "may still affect global domain "
                "orientation or protein dynamics."
            ),
            (
                "This analysis does not replace "
                "clash, bond-geometry, sequence-"
                "coverage, or confidence checks."
            ),
        ],
        "selected_pocket_ids": (
            selected_pocket_ids
        ),
        "missing_selected_pocket_ids": (
            missing_selected_pocket_ids
        ),
        "global_ramachandran_summary": (
            global_summary
        ),
        "global_goals_met": (
            global_goals_met
        ),
        "outlier_count": len(outliers),
        "unresolved_outlier_residues": (
            unresolved_outliers
        ),
        "inside_selected_box_outliers": (
            sorted(
                set(
                    inside_selected_box
                )
            )
        ),
        "near_selected_box_outliers": (
            sorted(
                set(
                    near_selected_box
                )
            )
        ),
        "selected_box_local_outliers": (
            selected_box_local_outliers
        ),
        "verdict": verdict,
        "outliers": outlier_results,
    }


def write_structure_pocket_quality(
    report: dict[str, Any],
    output_path: Path,
) -> Path:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path.write_text(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return output_path


def run_structure_pocket_quality(
    *,
    structure_path: Path,
    ramachandran_report_path: Path,
    pocket_definitions_path: Path,
    pocket_selection_summary_path: Path,
    output_dir: Path,
    near_threshold_angstrom: float = (
        DEFAULT_NEAR_BOX_THRESHOLD_ANGSTROM
    ),
) -> dict[str, Any]:
    report = evaluate_structure_pocket_quality(
        structure_path=structure_path,
        ramachandran_report_path=(
            ramachandran_report_path
        ),
        pocket_definitions_path=(
            pocket_definitions_path
        ),
        pocket_selection_summary_path=(
            pocket_selection_summary_path
        ),
        near_threshold_angstrom=(
            near_threshold_angstrom
        ),
    )

    output_path = (
        output_dir
        / "structure_pocket_quality.json"
    )

    write_structure_pocket_quality(
        report,
        output_path,
    )

    return {
        "report": report,
        "output_path": output_path,
    }
