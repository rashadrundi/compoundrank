from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PocketDefinition:
    mode: str
    center_x: float | None = None
    center_y: float | None = None
    center_z: float | None = None
    size_x: float | None = None
    size_y: float | None = None
    size_z: float | None = None
    autobox_ligand: Path | None = None
    source: str | None = None
    pocket_id: str = "pocket_01"
    pocket_rank: int = 1
    fpocket_score: float | None = None
    merged_from: tuple[str, ...] = field(
        default_factory=tuple
    )
    merge_distance: float | None = None

    def as_gnina_args(self, autobox_add: float = 4.0) -> list[str]:
        if self.mode == "autobox":
            if self.autobox_ligand is None:
                raise ValueError("Autobox mode requires autobox_ligand")
            return [
                "--autobox_ligand",
                str(self.autobox_ligand),
                "--autobox_add",
                str(autobox_add),
            ]

        if self.mode != "explicit":
            raise ValueError(f"Unsupported pocket mode: {self.mode}")

        values = (
            self.center_x,
            self.center_y,
            self.center_z,
            self.size_x,
            self.size_y,
            self.size_z,
        )
        if any(value is None for value in values):
            raise ValueError("Explicit pocket definition is incomplete")

        return [
            "--center_x",
            str(self.center_x),
            "--center_y",
            str(self.center_y),
            "--center_z",
            str(self.center_z),
            "--size_x",
            str(self.size_x),
            "--size_y",
            str(self.size_y),
            "--size_z",
            str(self.size_z),
        ]


@dataclass(frozen=True)
class PreparedReceptor:
    source_pdb: Path
    prepared_pdbqt: Path
    display_pdb: Path
    cache_key: str


@dataclass(frozen=True)
class PreparedLigand:
    name: str
    source_description: str
    source_sdf: Path
    prepared_pdbqt: Path
    cache_key: str


@dataclass
class PoseRecord:
    ligand_name: str
    seed: int
    pose_number: int
    molecule: Any
    cnn_score: float | None
    cnn_affinity: float | None
    minimized_affinity: float | None
    source_sdf: Path
    pocket_id: str = "pocket_01"
    pocket_rank: int = 1
    pocket_source: str | None = None
    fpocket_score: float | None = None


@dataclass
class PoseCluster:
    cluster_id: int
    representative: PoseRecord
    members: list[PoseRecord] = field(default_factory=list)
    seeds: set[int] = field(default_factory=set)
    member_count: int = 0
    valid_member_count: int = 0


@dataclass
class InteractionEvidence:
    contact_residues: list[str]
    polar_contact_candidates: list[str]
    hydrophobic_contact_residues: list[str]
    closest_residue_distance: float | None


@dataclass
class LigandResult:
    ligand: PreparedLigand
    clusters: list[PoseCluster]
    uncertainty: str
    uncertainty_reasons: list[str]
    top_score: float | None
