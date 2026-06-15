from __future__ import annotations

from pathlib import Path

from rdkit import Chem

from .chemistry import (
    choose_pose_to_source_mapping,
    load_first_sdf,
    load_sdf_records,
    parse_meeko_index_pairs,
    reconstruct_heavy_pose,
)
from .models import PocketDefinition, PoseRecord, PreparedLigand, PreparedReceptor
from .subprocess_utils import resolve_executable, run_command


def _float_property(molecule: Chem.Mol, name: str) -> float | None:
    if not molecule.HasProp(name):
        return None
    try:
        return float(molecule.GetProp(name))
    except ValueError:
        return None


def run_gnina_seed(
    receptor: PreparedReceptor,
    ligand: PreparedLigand,
    pocket: PocketDefinition,
    seed: int,
    work_dir: Path,
    *,
    exhaustiveness: int = 32,
    num_modes: int = 20,
    cnn_scoring: str = "refinement",
    autobox_add: float = 4.0,
    gnina_bin: str = "gnina",
    cpu: int | None = None,
    device: int | None = None,
) -> list[PoseRecord]:
    gnina = resolve_executable(gnina_bin, "GNINA")
    seed_dir = work_dir / ligand.name / f"seed_{seed}"
    seed_dir.mkdir(parents=True, exist_ok=True)
    output_sdf = seed_dir / "poses.sdf"

    command = [
        gnina,
        "--receptor",
        str(receptor.prepared_pdbqt),
        "--ligand",
        str(ligand.prepared_pdbqt),
        *pocket.as_gnina_args(autobox_add=autobox_add),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--seed",
        str(seed),
        "--cnn_scoring",
        cnn_scoring,
        "--out",
        str(output_sdf),
    ]
    if cpu is not None:
        command += ["--cpu", str(cpu)]
    if device is not None:
        command += ["--device", str(device)]

    completed = run_command(command)
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")

    if not output_sdf.is_file() or output_sdf.stat().st_size == 0:
        raise RuntimeError(f"GNINA did not create poses for {ligand.name}")

    raw_poses = load_sdf_records(output_sdf, sanitize=False)
    if not raw_poses:
        raise RuntimeError(f"GNINA pose SDF was unreadable: {output_sdf}")

    source = load_first_sdf(ligand.source_sdf, sanitize=True)
    pairs = parse_meeko_index_pairs(ligand.prepared_pdbqt)
    pose_to_source = choose_pose_to_source_mapping(pairs, source, raw_poses[0])

    records: list[PoseRecord] = []
    for pose_number, raw_pose in enumerate(raw_poses, start=1):
        reconstructed = reconstruct_heavy_pose(source, raw_pose, pose_to_source)
        reconstructed.SetProp("_Name", f"{ligand.name}_seed_{seed}_pose_{pose_number}")
        for property_name in ("CNNscore", "CNNaffinity", "minimizedAffinity"):
            value = _float_property(raw_pose, property_name)
            if value is not None:
                reconstructed.SetDoubleProp(property_name, value)
        reconstructed.SetIntProp("seed", seed)
        reconstructed.SetIntProp("pose_number", pose_number)
        records.append(
            PoseRecord(
                ligand_name=ligand.name,
                seed=seed,
                pose_number=pose_number,
                molecule=reconstructed,
                cnn_score=_float_property(raw_pose, "CNNscore"),
                cnn_affinity=_float_property(raw_pose, "CNNaffinity"),
                minimized_affinity=_float_property(raw_pose, "minimizedAffinity"),
                source_sdf=output_sdf,
            )
        )
    return records


def run_gnina_ensemble(
    receptor: PreparedReceptor,
    ligand: PreparedLigand,
    pocket: PocketDefinition,
    seeds: list[int],
    work_dir: Path,
    **kwargs: object,
) -> list[PoseRecord]:
    records: list[PoseRecord] = []
    for seed in seeds:
        print(f"\n[GNINA] {ligand.name}: seed {seed}")
        records.extend(
            run_gnina_seed(
                receptor,
                ligand,
                pocket,
                seed,
                work_dir,
                **kwargs,
            )
        )
    return records
