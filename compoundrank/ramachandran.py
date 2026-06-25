from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence

import gemmi


SCHEMA_VERSION = (
    "ramachandran_validation.v0.2"
)

CLASSIFICATION_METHOD = (
    "top8000_percentile_grid.v1"
)

REFERENCE_DATA_DIR = (
    Path(__file__).resolve().parent
    / "data"
    / "ramachandran"
)

FAVORED_CUTOFF = 0.02

ALLOWED_CUTOFFS = {
    "general": 0.0005,
    "glycine": 0.0010,
    "cis_proline": 0.0020,
    "trans_proline": 0.0010,
    "pre_proline": 0.0010,
    "isoleucine_valine": 0.0010,
}

GRID_FILES = {
    "general": (
        "rama8000-general-noGPIVpreP.data"
    ),
    "glycine": (
        "rama8000-gly-sym.data"
    ),
    "cis_proline": (
        "rama8000-cispro.data"
    ),
    "trans_proline": (
        "rama8000-transpro.data"
    ),
    "pre_proline": (
        "rama8000-prepro-noGP.data"
    ),
    "isoleucine_valine": (
        "rama8000-ileval-nopreP.data"
    ),
}

GRID_SHA256 = {
    "rama8000-general-noGPIVpreP.data": (
        "ccdbc6a201ca2510119e77b0dd169ca2"
        "cdea7f2c9ed348555a3b1129f8c2b00a"
    ),
    "rama8000-gly-sym.data": (
        "89c75a5ac036ff3309c51b46a2412f301"
        "37827abaff4606209ed41e04fcba637"
    ),
    "rama8000-cispro.data": (
        "143a3004668baefaf6749bd9c2acbe829d"
        "3b1a93059c74ba766b04f1e00f43ab"
    ),
    "rama8000-transpro.data": (
        "092b4c0bcd2fe846a063000c83d13cbe6"
        "2ff43b2dafd67ecb8a47042cce747ec"
    ),
    "rama8000-prepro-noGP.data": (
        "1cc10b92911d47ca775f0a8131a029b42"
        "634dbe26364cece4cfdc8c85b9c3fab"
    ),
    "rama8000-ileval-nopreP.data": (
        "567e1318128a50b8f44362427bf717d56"
        "4749936d22d5dbfd34e5b0c17819276"
    ),
}

NUMBER_PATTERN = re.compile(
    r"""
    [-+]?
    (?:
        (?:\d+\.\d*)
        |
        (?:\.\d+)
        |
        (?:\d+)
    )
    (?:[eE][-+]?\d+)?
    """,
    re.VERBOSE,
)


def _sha256(
    path: Path,
) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def _parse_grid(
    path: Path,
) -> dict[
    tuple[int, int],
    float,
]:
    grid: dict[
        tuple[int, int],
        float,
    ] = {}

    for line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        stripped = line.strip()

        if (
            not stripped
            or stripped.startswith("#")
        ):
            continue

        numbers = NUMBER_PATTERN.findall(
            stripped
        )

        if len(numbers) < 3:
            raise ValueError(
                "Malformed Top8000 row in "
                f"{path}: {line!r}"
            )

        phi_raw, psi_raw, score_raw = (
            numbers[-3:]
        )

        phi = int(
            round(
                float(phi_raw)
            )
        )

        psi = int(
            round(
                float(psi_raw)
            )
        )

        score = float(
            score_raw
        )

        grid[
            (
                phi,
                psi,
            )
        ] = score

    if not grid:
        raise ValueError(
            f"No grid points parsed: {path}"
        )

    return grid


@lru_cache(
    maxsize=4
)
def _load_grids_cached(
    directory_string: str,
    verify_checksums: bool,
) -> dict[
    str,
    dict[
        tuple[int, int],
        float,
    ],
]:
    directory = Path(
        directory_string
    )

    grids: dict[
        str,
        dict[
            tuple[int, int],
            float,
        ],
    ] = {}

    for category, filename in (
        GRID_FILES.items()
    ):
        path = directory / filename

        if not path.is_file():
            raise FileNotFoundError(
                path
            )

        if verify_checksums:
            actual = _sha256(
                path
            )

            expected = GRID_SHA256[
                filename
            ]

            if actual != expected:
                raise ValueError(
                    "Top8000 checksum "
                    f"mismatch for {path}. "
                    f"Expected {expected}; "
                    f"found {actual}."
                )

        grids[category] = (
            _parse_grid(
                path
            )
        )

    return grids


def load_top8000_grids(
    data_dir: Path | None = None,
    *,
    verify_checksums: bool = True,
) -> dict[
    str,
    dict[
        tuple[int, int],
        float,
    ],
]:
    directory = (
        Path(data_dir)
        if data_dir is not None
        else REFERENCE_DATA_DIR
    )

    return _load_grids_cached(
        str(
            directory.resolve()
        ),
        verify_checksums,
    )


def _normalize_angle(
    angle: float,
) -> float:
    return (
        (
            float(angle)
            + 180.0
        )
        % 360.0
        - 180.0
    )


def _grid_bin(
    angle: float,
) -> int:
    normalized = _normalize_angle(
        angle
    )

    index = math.floor(
        (
            normalized
            + 180.0
        )
        / 2.0
    )

    center = (
        -179
        + 2 * int(index)
    )

    return max(
        -179,
        min(
            179,
            center,
        ),
    )


def lookup_top8000_score(
    phi_degrees: float,
    psi_degrees: float,
    category: str,
    *,
    grids: (
        dict[
            str,
            dict[
                tuple[int, int],
                float,
            ],
        ]
        | None
    ) = None,
) -> dict[str, Any]:
    if category not in GRID_FILES:
        raise ValueError(
            "Unsupported Ramachandran "
            f"category: {category}"
        )

    effective_grids = (
        grids
        if grids is not None
        else load_top8000_grids()
    )

    phi_bin = _grid_bin(
        phi_degrees
    )

    psi_bin = _grid_bin(
        psi_degrees
    )

    score = effective_grids[
        category
    ].get(
        (
            phi_bin,
            psi_bin,
        ),
        0.0,
    )

    return {
        "score": float(
            score
        ),
        "phi_bin": phi_bin,
        "psi_bin": psi_bin,
        "grid_point_present": (
            (
                phi_bin,
                psi_bin,
            )
            in effective_grids[
                category
            ]
        ),
    }


def classify_score(
    category: str,
    score: float,
) -> str:
    if category not in (
        ALLOWED_CUTOFFS
    ):
        raise ValueError(
            "Unsupported Ramachandran "
            f"category: {category}"
        )

    value = float(
        score
    )

    if value >= FAVORED_CUTOFF:
        return "favored"

    if value >= ALLOWED_CUTOFFS[
        category
    ]:
        return "allowed"

    return "outlier"


def classify_ramachandran(
    phi_degrees: float,
    psi_degrees: float,
    category: str,
    *,
    grids: (
        dict[
            str,
            dict[
                tuple[int, int],
                float,
            ],
        ]
        | None
    ) = None,
) -> dict[str, Any]:
    lookup = lookup_top8000_score(
        phi_degrees,
        psi_degrees,
        category,
        grids=grids,
    )

    score = float(
        lookup["score"]
    )

    return {
        "classification": (
            classify_score(
                category,
                score,
            )
        ),
        "score": score,
        "score_percent": (
            score * 100.0
        ),
        "favored_cutoff": (
            FAVORED_CUTOFF
        ),
        "allowed_cutoff": (
            ALLOWED_CUTOFFS[
                category
            ]
        ),
        **lookup,
    }


def _residue_number(
    residue: gemmi.Residue,
) -> str:
    number = str(
        residue.seqid.num
    )

    insertion_code = str(
        residue.seqid.icode
    ).strip()

    if insertion_code:
        number += insertion_code

    return number


def _preferred_atom_position(
    residue: gemmi.Residue,
    atom_name: str,
) -> (
    tuple[
        float,
        float,
        float,
    ]
    | None
):
    candidates = [
        atom
        for atom in residue
        if (
            atom.name.strip().upper()
            == atom_name.upper()
        )
    ]

    if not candidates:
        return None

    def rank(
        atom: gemmi.Atom,
    ) -> tuple[int, float]:
        altloc = str(
            atom.altloc
        ).strip()

        if altloc in {
            "",
            "\x00",
        }:
            altloc_rank = 0
        elif altloc == "A":
            altloc_rank = 1
        else:
            altloc_rank = 2

        return (
            altloc_rank,
            -float(
                atom.occ
            ),
        )

    selected = min(
        candidates,
        key=rank,
    )

    return (
        float(
            selected.pos.x
        ),
        float(
            selected.pos.y
        ),
        float(
            selected.pos.z
        ),
    )


def _subtract(
    first: tuple[
        float,
        float,
        float,
    ],
    second: tuple[
        float,
        float,
        float,
    ],
) -> tuple[
    float,
    float,
    float,
]:
    return (
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    )


def _dot(
    first: tuple[
        float,
        float,
        float,
    ],
    second: tuple[
        float,
        float,
        float,
    ],
) -> float:
    return (
        first[0] * second[0]
        + first[1] * second[1]
        + first[2] * second[2]
    )


def _cross(
    first: tuple[
        float,
        float,
        float,
    ],
    second: tuple[
        float,
        float,
        float,
    ],
) -> tuple[
    float,
    float,
    float,
]:
    return (
        (
            first[1] * second[2]
            - first[2] * second[1]
        ),
        (
            first[2] * second[0]
            - first[0] * second[2]
        ),
        (
            first[0] * second[1]
            - first[1] * second[0]
        ),
    )


def _norm(
    vector: tuple[
        float,
        float,
        float,
    ],
) -> float:
    return math.sqrt(
        _dot(
            vector,
            vector,
        )
    )


def _scale(
    vector: tuple[
        float,
        float,
        float,
    ],
    factor: float,
) -> tuple[
    float,
    float,
    float,
]:
    return (
        vector[0] * factor,
        vector[1] * factor,
        vector[2] * factor,
    )


def _dihedral_degrees(
    first: tuple[
        float,
        float,
        float,
    ],
    second: tuple[
        float,
        float,
        float,
    ],
    third: tuple[
        float,
        float,
        float,
    ],
    fourth: tuple[
        float,
        float,
        float,
    ],
) -> float | None:
    first_bond = _subtract(
        first,
        second,
    )

    middle_bond = _subtract(
        third,
        second,
    )

    final_bond = _subtract(
        fourth,
        third,
    )

    middle_norm = _norm(
        middle_bond
    )

    if middle_norm == 0.0:
        return None

    middle_unit = _scale(
        middle_bond,
        1.0 / middle_norm,
    )

    first_projection = _subtract(
        first_bond,
        _scale(
            middle_unit,
            _dot(
                first_bond,
                middle_unit,
            ),
        ),
    )

    final_projection = _subtract(
        final_bond,
        _scale(
            middle_unit,
            _dot(
                final_bond,
                middle_unit,
            ),
        ),
    )

    if (
        _norm(
            first_projection
        )
        == 0.0
        or _norm(
            final_projection
        )
        == 0.0
    ):
        return None

    x_value = _dot(
        first_projection,
        final_projection,
    )

    y_value = _dot(
        _cross(
            middle_unit,
            first_projection,
        ),
        final_projection,
    )

    return math.degrees(
        math.atan2(
            y_value,
            x_value,
        )
    )


def _omega_degrees(
    previous_residue: gemmi.Residue,
    residue: gemmi.Residue,
) -> float | None:
    positions = (
        _preferred_atom_position(
            previous_residue,
            "CA",
        ),
        _preferred_atom_position(
            previous_residue,
            "C",
        ),
        _preferred_atom_position(
            residue,
            "N",
        ),
        _preferred_atom_position(
            residue,
            "CA",
        ),
    )

    if any(
        position is None
        for position in positions
    ):
        return None

    first, second, third, fourth = (
        positions
    )

    assert first is not None
    assert second is not None
    assert third is not None
    assert fourth is not None

    return _dihedral_degrees(
        first,
        second,
        third,
        fourth,
    )


def _residue_category(
    previous_residue: gemmi.Residue,
    residue: gemmi.Residue,
    next_residue: gemmi.Residue,
) -> tuple[
    str | None,
    float | None,
]:
    residue_name = (
        residue.name.upper()
    )

    next_name = (
        next_residue.name.upper()
    )

    if residue_name == "GLY":
        return (
            "glycine",
            None,
        )

    if residue_name == "PRO":
        omega = _omega_degrees(
            previous_residue,
            residue,
        )

        if omega is None:
            return (
                None,
                None,
            )

        normalized_omega = (
            _normalize_angle(
                omega
            )
        )

        if (
            -90.0
            < normalized_omega
            < 90.0
        ):
            return (
                "cis_proline",
                normalized_omega,
            )

        return (
            "trans_proline",
            normalized_omega,
        )

    if next_name == "PRO":
        return (
            "pre_proline",
            None,
        )

    if residue_name in {
        "ILE",
        "VAL",
    }:
        return (
            "isoleucine_valine",
            None,
        )

    return (
        "general",
        None,
    )


def _screening_flag(
    favored_fraction: float,
    outlier_fraction: float,
) -> str:
    favored_goal_met = (
        favored_fraction > 0.98
    )

    outlier_goal_met = (
        outlier_fraction < 0.002
    )

    if (
        favored_goal_met
        and outlier_goal_met
    ):
        return "meets_ramalyze_goals"

    if outlier_fraction >= 0.02:
        return "high_outlier_fraction"

    if outlier_fraction >= 0.002:
        return "elevated_outlier_fraction"

    return "favored_fraction_below_goal"


def validate_ramachandran(
    structure_path: Path,
    *,
    model_index: int = 0,
    chain_id: str | None = None,
) -> dict[str, Any]:
    source = Path(
        structure_path
    )

    if not source.is_file():
        raise FileNotFoundError(
            source
        )

    grids = load_top8000_grids()

    structure = gemmi.read_structure(
        str(source)
    )

    structure.setup_entities()

    if len(structure) == 0:
        raise ValueError(
            "Structure contains no models."
        )

    if not (
        0
        <= model_index
        < len(structure)
    ):
        raise ValueError(
            "Model index is outside the "
            "available model range."
        )

    model = structure[
        model_index
    ]

    model_identifier = getattr(
        model,
        "name",
        None,
    )

    if model_identifier is None:
        model_identifier = getattr(
            model,
            "num",
            None,
        )

    if callable(
        model_identifier
    ):
        model_identifier = (
            model_identifier()
        )

    if model_identifier is None:
        model_identifier = (
            model_index + 1
        )

    model_label = str(
        model_identifier
    )

    requested_chain = (
        str(chain_id).strip()
        if chain_id is not None
        else None
    )

    records: list[
        dict[str, Any]
    ] = []

    skipped: list[
        dict[str, Any]
    ] = []

    chain_counts: dict[
        str,
        Counter[str]
    ] = {}

    category_counts: dict[
        str,
        Counter[str]
    ] = {
        category: Counter()
        for category in GRID_FILES
    }

    total_polymer_residues = 0
    matched_chain = False

    for chain in model:
        current_chain = (
            chain.name.strip()
            or "_"
        )

        if (
            requested_chain is not None
            and current_chain
            != requested_chain
        ):
            continue

        matched_chain = True

        polymer = chain.get_polymer()

        total_polymer_residues += (
            len(polymer)
        )

        chain_counter: Counter[str] = (
            Counter()
        )

        chain_counts[
            current_chain
        ] = chain_counter

        for residue in polymer:
            residue_number = (
                _residue_number(
                    residue
                )
            )

            residue_id = (
                f"{residue.name}:"
                f"{current_chain}:"
                f"{residue_number}"
            )

            previous_residue = (
                chain.previous_residue(
                    residue
                )
            )

            next_residue = (
                chain.next_residue(
                    residue
                )
            )

            if (
                previous_residue is None
                or next_residue is None
            ):
                skipped.append(
                    {
                        "residue": (
                            residue_id
                        ),
                        "reason": (
                            "terminus_or_"
                            "chain_break"
                        ),
                    }
                )

                chain_counter[
                    "skipped"
                ] += 1

                continue

            category, omega_degrees = (
                _residue_category(
                    previous_residue,
                    residue,
                    next_residue,
                )
            )

            if category is None:
                skipped.append(
                    {
                        "residue": (
                            residue_id
                        ),
                        "reason": (
                            "proline_omega_"
                            "unavailable"
                        ),
                    }
                )

                chain_counter[
                    "skipped"
                ] += 1

                continue

            try:
                phi_radians, psi_radians = (
                    gemmi.calculate_phi_psi(
                        previous_residue,
                        residue,
                        next_residue,
                    )
                )
            except Exception as error:
                skipped.append(
                    {
                        "residue": (
                            residue_id
                        ),
                        "reason": (
                            "torsion_"
                            "calculation_failed"
                        ),
                        "error": str(
                            error
                        ),
                    }
                )

                chain_counter[
                    "skipped"
                ] += 1

                continue

            phi_degrees = math.degrees(
                phi_radians
            )

            psi_degrees = math.degrees(
                psi_radians
            )

            if (
                not math.isfinite(
                    phi_degrees
                )
                or not math.isfinite(
                    psi_degrees
                )
            ):
                skipped.append(
                    {
                        "residue": (
                            residue_id
                        ),
                        "reason": (
                            "missing_backbone_"
                            "coordinates"
                        ),
                    }
                )

                chain_counter[
                    "skipped"
                ] += 1

                continue

            classification = (
                classify_ramachandran(
                    phi_degrees,
                    psi_degrees,
                    category,
                    grids=grids,
                )
            )

            result_class = str(
                classification[
                    "classification"
                ]
            )

            chain_counter[
                result_class
            ] += 1

            chain_counter[
                "evaluable"
            ] += 1

            category_counts[
                category
            ][
                result_class
            ] += 1

            category_counts[
                category
            ][
                "evaluable"
            ] += 1

            records.append(
                {
                    "model": model_label,
                    "chain": (
                        current_chain
                    ),
                    "residue_name": (
                        residue.name
                    ),
                    "residue_number": (
                        residue_number
                    ),
                    "residue": (
                        residue_id
                    ),
                    "category": (
                        category
                    ),
                    "omega_degrees": (
                        omega_degrees
                    ),
                    "phi_degrees": (
                        phi_degrees
                    ),
                    "psi_degrees": (
                        psi_degrees
                    ),
                    **classification,
                }
            )

    if (
        requested_chain is not None
        and not matched_chain
    ):
        raise ValueError(
            "Requested chain was not "
            f"found: {requested_chain}"
        )

    classification_counts = Counter(
        str(
            record[
                "classification"
            ]
        )
        for record in records
    )

    evaluable_count = len(
        records
    )

    favored_count = int(
        classification_counts[
            "favored"
        ]
    )

    allowed_count = int(
        classification_counts[
            "allowed"
        ]
    )

    outlier_count = int(
        classification_counts[
            "outlier"
        ]
    )

    def fraction(
        count: int,
    ) -> float:
        if evaluable_count == 0:
            return 0.0

        return (
            count
            / evaluable_count
        )

    favored_fraction = fraction(
        favored_count
    )

    allowed_fraction = fraction(
        allowed_count
    )

    outlier_fraction = fraction(
        outlier_count
    )

    per_chain: dict[
        str,
        dict[str, Any],
    ] = {}

    for current_chain, counts in sorted(
        chain_counts.items()
    ):
        chain_evaluable = int(
            counts["evaluable"]
        )

        chain_outliers = int(
            counts["outlier"]
        )

        per_chain[
            current_chain
        ] = {
            "evaluable_residues": (
                chain_evaluable
            ),
            "favored": int(
                counts["favored"]
            ),
            "allowed": int(
                counts["allowed"]
            ),
            "outliers": (
                chain_outliers
            ),
            "skipped": int(
                counts["skipped"]
            ),
            "outlier_fraction": (
                chain_outliers
                / chain_evaluable
                if chain_evaluable
                else 0.0
            ),
        }

    per_category: dict[
        str,
        dict[str, Any],
    ] = {}

    for category, counts in (
        category_counts.items()
    ):
        category_evaluable = int(
            counts["evaluable"]
        )

        per_category[category] = {
            "evaluable_residues": (
                category_evaluable
            ),
            "favored": int(
                counts["favored"]
            ),
            "allowed": int(
                counts["allowed"]
            ),
            "outliers": int(
                counts["outlier"]
            ),
            "allowed_cutoff": (
                ALLOWED_CUTOFFS[
                    category
                ]
            ),
        }

    return {
        "schema_version": (
            SCHEMA_VERSION
        ),
        "status": (
            "complete"
            if evaluable_count > 0
            else (
                "insufficient_geometry"
            )
        ),
        "selection_mode": (
            "report_only"
        ),
        "classification_method": (
            CLASSIFICATION_METHOD
        ),
        "reference_data": {
            "dataset": (
                "Top8000 Ramachandran "
                "percentile contour grids"
            ),
            "origin": (
                "Richardson Laboratory "
                "reference_data repository"
            ),
            "license": (
                "CC BY 4.0"
            ),
            "directory": str(
                REFERENCE_DATA_DIR
            ),
            "favored_cutoff": (
                FAVORED_CUTOFF
            ),
            "allowed_cutoffs": (
                ALLOWED_CUTOFFS
            ),
            "grid_files": (
                GRID_FILES
            ),
            "grid_sha256": (
                GRID_SHA256
            ),
        },
        "source_structure": str(
            source.resolve()
        ),
        "structure_name": (
            structure.name
        ),
        "model_index": (
            model_index
        ),
        "model_name": (
            model_label
        ),
        "requested_chain": (
            requested_chain
        ),
        "total_polymer_residues": (
            total_polymer_residues
        ),
        "evaluable_residues": (
            evaluable_count
        ),
        "skipped_residues": len(
            skipped
        ),
        "summary": {
            "favored": (
                favored_count
            ),
            "allowed": (
                allowed_count
            ),
            "outliers": (
                outlier_count
            ),
            "favored_fraction": (
                favored_fraction
            ),
            "allowed_fraction": (
                allowed_fraction
            ),
            "outlier_fraction": (
                outlier_fraction
            ),
            "favored_goal": (
                "> 98%"
            ),
            "outlier_goal": (
                "< 0.2%"
            ),
            "favored_goal_met": (
                favored_fraction > 0.98
            ),
            "outlier_goal_met": (
                outlier_fraction < 0.002
            ),
            "screening_flag": (
                _screening_flag(
                    favored_fraction,
                    outlier_fraction,
                )
                if evaluable_count
                else (
                    "insufficient_geometry"
                )
            ),
        },
        "per_chain": (
            per_chain
        ),
        "per_category": (
            per_category
        ),
        "residues": records,
        "skipped": skipped,
        "limitations": [
            (
                "This release is "
                "report-only and does not "
                "reject or rerank structures."
            ),
            (
                "Grid values are looked up "
                "at the corresponding 2-degree "
                "Top8000 bin; sparse missing "
                "bins have score zero."
            ),
            (
                "Residues adjacent to chain "
                "breaks or lacking complete "
                "backbone coordinates are "
                "not evaluated."
            ),
            (
                "Ramachandran validation is "
                "one component of structure "
                "quality and does not replace "
                "clash, rotamer, bond, angle, "
                "density, or dynamics checks."
            ),
        ],
    }


def write_ramachandran_outputs(
    report: dict[str, Any],
    output_dir: Path,
) -> dict[str, str]:
    destination = Path(
        output_dir
    )

    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    json_path = (
        destination
        / "ramachandran_validation.json"
    )

    csv_path = (
        destination
        / "ramachandran_residues.csv"
    )

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    fieldnames = [
        "model",
        "chain",
        "residue_name",
        "residue_number",
        "residue",
        "category",
        "omega_degrees",
        "phi_degrees",
        "psi_degrees",
        "phi_bin",
        "psi_bin",
        "grid_point_present",
        "score",
        "score_percent",
        "classification",
        "favored_cutoff",
        "allowed_cutoff",
    ]

    with csv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for record in report[
            "residues"
        ]:
            writer.writerow(
                {
                    key: record.get(
                        key
                    )
                    for key in fieldnames
                }
            )

    return {
        "json": str(
            json_path
        ),
        "csv": str(
            csv_path
        ),
    }


def build_cli_parser() -> (
    argparse.ArgumentParser
):
    parser = argparse.ArgumentParser(
        description=(
            "Calculate backbone phi/psi "
            "angles and classify them with "
            "the Top8000 Ramachandran "
            "percentile contour grids."
        )
    )

    parser.add_argument(
        "--input-pdb",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--model-index",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--chain",
        default=None,
    )

    return parser


def main(
    argv: Sequence[str] | None = None,
) -> int:
    parser = build_cli_parser()

    arguments = parser.parse_args(
        argv
    )

    report = validate_ramachandran(
        arguments.input_pdb,
        model_index=(
            arguments.model_index
        ),
        chain_id=arguments.chain,
    )

    outputs = (
        write_ramachandran_outputs(
            report,
            arguments.output_dir,
        )
    )

    print(
        json.dumps(
            {
                "status": (
                    report["status"]
                ),
                "selection_mode": (
                    report[
                        "selection_mode"
                    ]
                ),
                "classification_method": (
                    report[
                        "classification_method"
                    ]
                ),
                "source_structure": (
                    report[
                        "source_structure"
                    ]
                ),
                "evaluable_residues": (
                    report[
                        "evaluable_residues"
                    ]
                ),
                "summary": (
                    report["summary"]
                ),
                "outputs": outputs,
            },
            indent=2,
            sort_keys=True,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
