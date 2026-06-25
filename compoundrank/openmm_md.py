from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

from .structure_ensemble import (
    build_structure_ensemble,
)


SCHEMA_VERSION = "openmm_md_ensemble.v0.1"


@dataclass(frozen=True)
class OpenMMMDConfig:
    temperature_kelvin: float = 300.0
    friction_per_ps: float = 1.0
    timestep_fs: float = 2.0
    equilibration_steps: int = 5000
    production_steps: int = 25000
    snapshot_interval: int = 5000
    seed: int = 20260625
    ph: float = 7.4
    hydrogen_mass_amu: float = 3.0
    minimization_tolerance_kj_mol_nm: float = 10.0
    minimization_max_iterations: int = 2000
    forcefield_file: str = "amber14-all.xml"
    implicit_solvent_file: str = "implicit/obc2.xml"


def validate_config(
    config: OpenMMMDConfig,
) -> None:
    positive_floats = {
        "temperature_kelvin": config.temperature_kelvin,
        "friction_per_ps": config.friction_per_ps,
        "timestep_fs": config.timestep_fs,
        "ph": config.ph,
        "hydrogen_mass_amu": config.hydrogen_mass_amu,
        "minimization_tolerance_kj_mol_nm": (
            config.minimization_tolerance_kj_mol_nm
        ),
    }

    for name, value in positive_floats.items():
        if (
            not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError(
                f"{name} must be a finite "
                "positive value"
            )

    positive_integers = {
        "equilibration_steps": (
            config.equilibration_steps
        ),
        "production_steps": (
            config.production_steps
        ),
        "snapshot_interval": (
            config.snapshot_interval
        ),
        "minimization_max_iterations": (
            config.minimization_max_iterations
        ),
    }

    for name, value in positive_integers.items():
        if value <= 0:
            raise ValueError(
                f"{name} must be greater "
                "than zero"
            )

    if (
        config.production_steps
        % config.snapshot_interval
        != 0
    ):
        raise ValueError(
            "production_steps must be evenly "
            "divisible by snapshot_interval"
        )


def choose_platform_name(
    available_platforms: Sequence[str],
    requested_platform: str,
) -> str:
    available = list(
        available_platforms
    )

    if requested_platform != "auto":
        if requested_platform not in available:
            raise RuntimeError(
                "Requested OpenMM platform "
                f"is unavailable: "
                f"{requested_platform}. "
                f"Available platforms: "
                f"{available}"
            )

        return requested_platform

    for candidate in (
        "CUDA",
        "OpenCL",
        "CPU",
        "Reference",
    ):
        if candidate in available:
            return candidate

    raise RuntimeError(
        "No OpenMM computation platform "
        "is available"
    )


def build_pdb2pqr_command(
    *,
    pdb2pqr_bin: str,
    receptor_pdb: Path,
    output_pdb: Path,
    output_pqr: Path,
    ph: float,
) -> list[str]:
    return [
        pdb2pqr_bin,
        "--ff=AMBER",
        "--keep-chain",
        "--titration-state-method=propka",
        f"--with-ph={ph}",
        "--pdb-output",
        str(output_pdb),
        str(receptor_pdb),
        str(output_pqr),
    ]


def _load_openmm() -> tuple[
    Any,
    Any,
    Any,
]:
    try:
        import openmm as mm
        from openmm import app, unit

    except ImportError as error:
        raise RuntimeError(
            "OpenMM is not installed in the "
            "active Python environment"
        ) from error

    return mm, app, unit


def _run_pdb2pqr(
    *,
    pdb2pqr_bin: str,
    receptor_pdb: Path,
    output_pdb: Path,
    output_pqr: Path,
    log_path: Path,
    ph: float,
) -> None:
    command = build_pdb2pqr_command(
        pdb2pqr_bin=pdb2pqr_bin,
        receptor_pdb=receptor_pdb,
        output_pdb=output_pdb,
        output_pqr=output_pqr,
        ph=ph,
    )

    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    log_path.write_text(
        result.stdout,
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RuntimeError(
            "PDB2PQR failed with exit code "
            f"{result.returncode}. "
            f"See {log_path}"
        )

    for output in (
        output_pdb,
        output_pqr,
    ):
        if (
            not output.is_file()
            or output.stat().st_size == 0
        ):
            raise RuntimeError(
                "PDB2PQR did not create a "
                f"usable output: {output}"
            )


def _potential_energy_kj_mol(
    state: Any,
    unit: Any,
) -> float:
    return float(
        state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )
    )


def _write_pdb(
    *,
    app: Any,
    topology: Any,
    positions: Any,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        app.PDBFile.writeFile(
            topology,
            positions,
            handle,
            keepIds=True,
        )


def _write_json_report(
    report: dict[str, Any],
    output_dir: Path,
) -> Path:
    report_path = (
        output_dir
        / "openmm_md_report.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return report_path


def run_openmm_md_ensemble(
    *,
    receptor_pdb: Path,
    output_dir: Path,
    config: OpenMMMDConfig,
    platform_name: str = "auto",
    device_index: str | None = None,
    pdb2pqr_bin: str = "pdb2pqr",
    overwrite: bool = False,
) -> dict[str, Any]:
    validate_config(
        config
    )

    receptor = (
        Path(receptor_pdb)
        .expanduser()
        .resolve()
    )

    destination = (
        Path(output_dir)
        .expanduser()
        .resolve()
    )

    if (
        not receptor.is_file()
        or receptor.stat().st_size == 0
    ):
        raise FileNotFoundError(
            receptor
        )

    if (
        destination.exists()
        and any(destination.iterdir())
        and not overwrite
    ):
        raise FileExistsError(
            "OpenMM output directory is "
            f"not empty: {destination}"
        )

    if (
        overwrite
        and destination.exists()
    ):
        shutil.rmtree(
            destination
        )

    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    prepared_dir = (
        destination / "prepared"
    )

    snapshot_dir = (
        destination
        / "production_snapshots"
    )

    prepared_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    snapshot_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    prepared_pdb = (
        prepared_dir
        / "receptor_pdb2pqr.pdb"
    )

    prepared_pqr = (
        prepared_dir
        / "receptor_pdb2pqr.pqr"
    )

    pdb2pqr_log = (
        prepared_dir
        / "pdb2pqr.log"
    )

    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "started",
        "selection_mode": "report_only",
        "source_receptor": str(
            receptor
        ),
        "output_dir": str(
            destination
        ),
        "requested_platform": (
            platform_name
        ),
        "config": asdict(
            config
        ),
    }

    try:
        _run_pdb2pqr(
            pdb2pqr_bin=pdb2pqr_bin,
            receptor_pdb=receptor,
            output_pdb=prepared_pdb,
            output_pqr=prepared_pqr,
            log_path=pdb2pqr_log,
            ph=config.ph,
        )

        mm, app, unit = (
            _load_openmm()
        )

        available_platforms = [
            mm.Platform
            .getPlatform(index)
            .getName()
            for index in range(
                mm.Platform.getNumPlatforms()
            )
        ]

        selected_platform_name = (
            choose_platform_name(
                available_platforms,
                platform_name,
            )
        )

        platform = (
            mm.Platform.getPlatformByName(
                selected_platform_name
            )
        )

        property_names = set(
            platform.getPropertyNames()
        )

        platform_properties: dict[
            str,
            str,
        ] = {}

        if (
            selected_platform_name == "CUDA"
            and "Precision" in property_names
        ):
            platform_properties[
                "Precision"
            ] = "mixed"

        if (
            device_index is not None
            and "DeviceIndex" in property_names
        ):
            platform_properties[
                "DeviceIndex"
            ] = device_index

        pdb = app.PDBFile(
            str(prepared_pdb)
        )

        forcefield = app.ForceField(
            config.forcefield_file,
            config.implicit_solvent_file,
        )

        unmatched = (
            forcefield.getUnmatchedResidues(
                pdb.topology
            )
        )

        unmatched_rows = [
            {
                "index": residue.index,
                "name": residue.name,
                "id": residue.id,
                "chain": residue.chain.id,
            }
            for residue in unmatched
        ]

        if unmatched_rows:
            raise RuntimeError(
                "OpenMM force-field templates "
                "remain unmatched: "
                f"{unmatched_rows}"
            )

        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
            hydrogenMass=(
                config.hydrogen_mass_amu
                * unit.amu
            ),
        )

        integrator = (
            mm.LangevinMiddleIntegrator(
                (
                    config.temperature_kelvin
                    * unit.kelvin
                ),
                (
                    config.friction_per_ps
                    / unit.picosecond
                ),
                (
                    config.timestep_fs
                    * unit.femtoseconds
                ),
            )
        )

        integrator.setRandomNumberSeed(
            config.seed
        )

        simulation = app.Simulation(
            pdb.topology,
            system,
            integrator,
            platform,
            platform_properties,
        )

        simulation.context.setPositions(
            pdb.positions
        )

        initial_state = (
            simulation.context.getState(
                getEnergy=True
            )
        )

        initial_energy = (
            _potential_energy_kj_mol(
                initial_state,
                unit,
            )
        )

        simulation.minimizeEnergy(
            tolerance=(
                config
                .minimization_tolerance_kj_mol_nm
                * unit.kilojoule_per_mole
                / unit.nanometer
            ),
            maxIterations=(
                config
                .minimization_max_iterations
            ),
        )

        minimized_state = (
            simulation.context.getState(
                getEnergy=True,
                getPositions=True,
            )
        )

        minimized_energy = (
            _potential_energy_kj_mol(
                minimized_state,
                unit,
            )
        )

        minimized_pdb = (
            destination
            / "receptor_minimized.pdb"
        )

        _write_pdb(
            app=app,
            topology=pdb.topology,
            positions=(
                minimized_state
                .getPositions()
            ),
            output_path=minimized_pdb,
        )

        simulation.context.setVelocitiesToTemperature(
            (
                config.temperature_kelvin
                * unit.kelvin
            ),
            config.seed,
        )

        simulation.step(
            config.equilibration_steps
        )

        equilibrated_state = (
            simulation.context.getState(
                getEnergy=True,
                getPositions=True,
            )
        )

        equilibrated_energy = (
            _potential_energy_kj_mol(
                equilibrated_state,
                unit,
            )
        )

        equilibrated_pdb = (
            destination
            / "receptor_equilibrated.pdb"
        )

        _write_pdb(
            app=app,
            topology=pdb.topology,
            positions=(
                equilibrated_state
                .getPositions()
            ),
            output_path=equilibrated_pdb,
        )

        state_log = (
            destination
            / "production_state.csv"
        )

        trajectory_path = (
            destination
            / "production_trajectory.dcd"
        )

        simulation.reporters.append(
            app.StateDataReporter(
                str(state_log),
                config.snapshot_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                separator=",",
            )
        )

        simulation.reporters.append(
            app.DCDReporter(
                str(trajectory_path),
                config.snapshot_interval,
            )
        )

        snapshot_count = (
            config.production_steps
            // config.snapshot_interval
        )

        snapshot_paths: list[
            Path
        ] = []

        snapshot_rows: list[
            dict[str, Any]
        ] = []

        for snapshot_index in range(
            1,
            snapshot_count + 1,
        ):
            simulation.step(
                config.snapshot_interval
            )

            state = (
                simulation.context.getState(
                    getEnergy=True,
                    getPositions=True,
                )
            )

            production_step = (
                snapshot_index
                * config.snapshot_interval
            )

            total_step = (
                config.equilibration_steps
                + production_step
            )

            production_time_ps = (
                production_step
                * config.timestep_fs
                / 1000.0
            )

            total_time_ps = (
                total_step
                * config.timestep_fs
                / 1000.0
            )

            snapshot_path = (
                snapshot_dir
                / (
                    f"snapshot_"
                    f"{snapshot_index:04d}.pdb"
                )
            )

            _write_pdb(
                app=app,
                topology=pdb.topology,
                positions=(
                    state.getPositions()
                ),
                output_path=snapshot_path,
            )

            snapshot_paths.append(
                snapshot_path
            )

            snapshot_rows.append(
                {
                    "snapshot_index": (
                        snapshot_index
                    ),
                    "production_step": (
                        production_step
                    ),
                    "total_step": total_step,
                    "production_time_ps": (
                        production_time_ps
                    ),
                    "total_time_ps": (
                        total_time_ps
                    ),
                    (
                        "potential_energy_"
                        "kj_mol"
                    ): (
                        _potential_energy_kj_mol(
                            state,
                            unit,
                        )
                    ),
                    "pdb": str(
                        snapshot_path
                    ),
                }
            )

        checkpoint_path = (
            destination
            / "final_checkpoint.chk"
        )

        simulation.saveCheckpoint(
            str(checkpoint_path)
        )

        ensemble_dir = (
            destination
            / "structure_ensemble"
        )

        ensemble_manifest = (
            build_structure_ensemble(
                reference_pdb=prepared_pdb,
                snapshot_pdbs=(
                    snapshot_paths
                ),
                output_dir=ensemble_dir,
                source_engine="openmm",
                overwrite=True,
            )
        )

        report.update(
            {
                "status": "complete",
                "openmm_version": (
                    mm.version.version
                ),
                "available_platforms": (
                    available_platforms
                ),
                "platform": (
                    selected_platform_name
                ),
                "platform_properties": (
                    platform_properties
                ),
                "prepared_receptor": str(
                    prepared_pdb
                ),
                "prepared_pqr": str(
                    prepared_pqr
                ),
                "pdb2pqr_log": str(
                    pdb2pqr_log
                ),
                "residue_count": (
                    pdb.topology
                    .getNumResidues()
                ),
                "atom_count": (
                    pdb.topology
                    .getNumAtoms()
                ),
                "unmatched_residue_count": 0,
                (
                    "initial_potential_"
                    "energy_kj_mol"
                ): initial_energy,
                (
                    "minimized_potential_"
                    "energy_kj_mol"
                ): minimized_energy,
                (
                    "equilibrated_potential_"
                    "energy_kj_mol"
                ): equilibrated_energy,
                "equilibration_time_ps": (
                    config.equilibration_steps
                    * config.timestep_fs
                    / 1000.0
                ),
                "production_time_ps": (
                    config.production_steps
                    * config.timestep_fs
                    / 1000.0
                ),
                "snapshot_count": (
                    snapshot_count
                ),
                "snapshots": snapshot_rows,
                "artifacts": {
                    "minimized_pdb": str(
                        minimized_pdb
                    ),
                    "equilibrated_pdb": str(
                        equilibrated_pdb
                    ),
                    "state_log": str(
                        state_log
                    ),
                    "trajectory_dcd": str(
                        trajectory_path
                    ),
                    "checkpoint": str(
                        checkpoint_path
                    ),
                    (
                        "structure_ensemble_"
                        "json"
                    ): (
                        ensemble_manifest[
                            "outputs"
                        ]["json"]
                    ),
                    (
                        "structure_ensemble_"
                        "csv"
                    ): (
                        ensemble_manifest[
                            "outputs"
                        ]["csv"]
                    ),
                },
                "limitations": [
                    (
                        "This initial backend "
                        "uses implicit OBC2 "
                        "solvent rather than an "
                        "explicit water box."
                    ),
                    (
                        "The ensemble is "
                        "report-only and does "
                        "not yet alter docking "
                        "or pocket selection."
                    ),
                    (
                        "Short validation runs "
                        "do not establish "
                        "biological convergence."
                    ),
                ],
            }
        )

        report_path = (
            _write_json_report(
                report,
                destination,
            )
        )

        report[
            "report_path"
        ] = str(
            report_path
        )

        _write_json_report(
            report,
            destination,
        )

        return report

    except Exception as error:
        report.update(
            {
                "status": "failed",
                "error_type": (
                    type(error).__name__
                ),
                "error": str(error),
                "traceback": (
                    traceback.format_exc()
                ),
            }
        )

        _write_json_report(
            report,
            destination,
        )

        raise


def build_cli_parser() -> (
    argparse.ArgumentParser
):
    parser = argparse.ArgumentParser(
        description=(
            "Generate an OpenMM receptor "
            "conformational ensemble and "
            "portable CompoundRank manifest."
        )
    )

    parser.add_argument(
        "--receptor",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--platform",
        choices=(
            "auto",
            "CUDA",
            "OpenCL",
            "CPU",
            "Reference",
        ),
        default="auto",
    )

    parser.add_argument(
        "--device-index",
        default=None,
    )

    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
    )

    parser.add_argument(
        "--friction",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--timestep-fs",
        type=float,
        default=2.0,
    )

    parser.add_argument(
        "--equilibration-steps",
        type=int,
        default=5000,
    )

    parser.add_argument(
        "--production-steps",
        type=int,
        default=25000,
    )

    parser.add_argument(
        "--snapshot-interval",
        type=int,
        default=5000,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20260625,
    )

    parser.add_argument(
        "--ph",
        type=float,
        default=7.4,
    )

    parser.add_argument(
        "--pdb2pqr-bin",
        default="pdb2pqr",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    return parser


def main() -> int:
    parser = build_cli_parser()
    args = parser.parse_args()

    config = OpenMMMDConfig(
        temperature_kelvin=(
            args.temperature
        ),
        friction_per_ps=(
            args.friction
        ),
        timestep_fs=(
            args.timestep_fs
        ),
        equilibration_steps=(
            args.equilibration_steps
        ),
        production_steps=(
            args.production_steps
        ),
        snapshot_interval=(
            args.snapshot_interval
        ),
        seed=args.seed,
        ph=args.ph,
    )

    report = run_openmm_md_ensemble(
        receptor_pdb=args.receptor,
        output_dir=args.output_dir,
        config=config,
        platform_name=args.platform,
        device_index=(
            args.device_index
        ),
        pdb2pqr_bin=(
            args.pdb2pqr_bin
        ),
        overwrite=args.overwrite,
    )

    print(
        json.dumps(
            {
                "status": (
                    report["status"]
                ),
                "platform": (
                    report["platform"]
                ),
                "snapshot_count": (
                    report[
                        "snapshot_count"
                    ]
                ),
                "report": (
                    report[
                        "report_path"
                    ]
                ),
                "ensemble_manifest": (
                    report[
                        "artifacts"
                    ][
                        "structure_ensemble_json"
                    ]
                ),
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
