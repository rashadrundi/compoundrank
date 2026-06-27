from __future__ import annotations

import shutil
from pathlib import Path

from .models import PreparedReceptor
from .paths import content_cache_key
from .subprocess_utils import resolve_executable, run_command


def _split_pdb_by_chain(source: Path, output_dir: Path) -> list[Path]:
    chain_lines: dict[str, list[str]] = {}
    for line in source.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21].strip() or "_"
        chain_lines.setdefault(chain, []).append(line)
    paths: list[Path] = []
    for index, (chain, lines) in enumerate(sorted(chain_lines.items()), start=1):
        path = output_dir / f"receptor_chain_{index:02d}_{chain}.pdb"
        path.write_text("\n".join([*lines, "TER", "END"]) + "\n", encoding="utf-8")
        paths.append(path)
    return paths


def _merge_pdbqt(inputs: list[Path], output: Path) -> None:
    lines: list[str] = []
    serial = 1
    for input_index, path in enumerate(inputs):
        for line in path.read_text(errors="replace").splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            lines.append(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        if input_index < len(inputs) - 1:
            lines.append("TER")
    lines.append("END")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run_meeko_receptor(
    meeko: str,
    input_pdb: Path,
    output_basename: Path,
    output_pdbqt: Path,
) -> None:
    run_command(
        [
            meeko,
            "--read_pdb",
            str(input_pdb),
            "--charge_model",
            "gasteiger",
            "-o",
            str(output_basename),
            "-p",
            str(output_pdbqt),
            "-j",
            str(output_basename.with_suffix(".json")),
        ]
    )



def _pdb_atom_element(
    line: str,
) -> str:
    """Return a PDB atom element conservatively."""

    if len(line) >= 78:
        explicit = (
            line[76:78]
            .strip()
            .upper()
        )

        if explicit:
            return explicit

    atom_name = (
        line[12:16].strip()
        if len(line) >= 16
        else ""
    )

    for character in atom_name:
        if character.isalpha():
            return character.upper()

    return ""


def _write_pdb2pqr_input(
    source_pdb: Path,
    destination_pdb: Path,
) -> int:
    """Write a heavy-atom PDB for PDB2PQR.

    MD/OpenMM outputs already contain explicit
    hydrogens. PDB2PQR must regenerate those
    hydrogens itself so protonation and debumping
    are internally consistent.
    """

    source = Path(source_pdb)
    destination = Path(
        destination_pdb
    )

    lines = source.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    retained: list[str] = []
    removed_hydrogens = 0
    retained_atom_count = 0

    for line in lines:
        if line.startswith(
            (
                "ATOM  ",
                "HETATM",
            )
        ):
            if _pdb_atom_element(
                line
            ) in {
                "H",
                "D",
            }:
                removed_hydrogens += 1
                continue

            retained_atom_count += 1

        retained.append(line)

    if retained_atom_count == 0:
        raise ValueError(
            "Receptor contains no heavy atoms "
            f"after hydrogen removal: {source}"
        )

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    destination.write_text(
        "\n".join(retained) + "\n",
        encoding="utf-8",
    )

    return removed_hydrogens


def prepare_receptor(
    source_pdb: Path,
    cache_root: Path,
    *,
    ph: float = 7.4,
    pdb2pqr_bin: str = "pdb2pqr",
    meeko_receptor_bin: str = "mk_prepare_receptor.py",
) -> PreparedReceptor:
    cache_key = content_cache_key(
        source_pdb,
        f"ph={ph}",
        "receptor-v4-heavy-input",
    )
    cache_dir = cache_root / "receptors" / cache_key
    prepared_pdbqt = cache_dir / "receptor_prepared.pdbqt"
    protonated_pdb = cache_dir / "receptor_protonated.pdb"
    pqr_path = cache_dir / "receptor.pqr"

    if prepared_pdbqt.is_file() and protonated_pdb.is_file():
        return PreparedReceptor(
            source_pdb=source_pdb,
            prepared_pdbqt=prepared_pdbqt,
            display_pdb=protonated_pdb,
            cache_key=cache_key,
        )

    cache_dir.mkdir(parents=True, exist_ok=True)
    copied_source = (
        cache_dir
        / "receptor_source.pdb"
    )

    shutil.copy2(
        source_pdb,
        copied_source,
    )

    pdb2pqr_input = (
        cache_dir
        / "receptor_pdb2pqr_input.pdb"
    )

    removed_hydrogens = (
        _write_pdb2pqr_input(
            copied_source,
            pdb2pqr_input,
        )
    )

    if removed_hydrogens:
        print(
            "[RECEPTOR PREPARATION] "
            f"Removed {removed_hydrogens} "
            "pre-existing hydrogen/deuterium "
            "atoms before PDB2PQR"
        )

    pdb2pqr = resolve_executable(pdb2pqr_bin, "PDB2PQR")
    meeko = resolve_executable(meeko_receptor_bin, "Meeko receptor preparation")

    run_command(
        [
            pdb2pqr,
            "--ff=AMBER",
            "--keep-chain",
            "--titration-state-method=propka",
            f"--with-ph={ph}",
            "--pdb-output",
            str(protonated_pdb),
            str(pdb2pqr_input),
            str(pqr_path),
        ]
    )

    try:
        _run_meeko_receptor(
            meeko,
            protonated_pdb,
            cache_dir / "receptor_prepared",
            prepared_pdbqt,
        )
    except RuntimeError as whole_error:
        # Some multimers are interpreted as having spurious inter-chain bonds.
        # Preparing chains separately and merging preserves chain coordinates and
        # avoids that parser failure. This is a fallback, not the first choice.
        chain_dir = cache_dir / "chains"
        chain_dir.mkdir(parents=True, exist_ok=True)
        chain_pdbs = _split_pdb_by_chain(protonated_pdb, chain_dir)
        if len(chain_pdbs) <= 1:
            raise whole_error
        prepared_chains: list[Path] = []
        for index, chain_pdb in enumerate(chain_pdbs, start=1):
            chain_pdbqt = chain_dir / f"chain_{index:02d}_prepared.pdbqt"
            _run_meeko_receptor(
                meeko,
                chain_pdb,
                chain_dir / f"chain_{index:02d}_prepared",
                chain_pdbqt,
            )
            prepared_chains.append(chain_pdbqt)
        _merge_pdbqt(prepared_chains, prepared_pdbqt)

    if not prepared_pdbqt.is_file() or prepared_pdbqt.stat().st_size == 0:
        raise RuntimeError("Meeko did not create a receptor PDBQT")

    return PreparedReceptor(
        source_pdb=source_pdb,
        prepared_pdbqt=prepared_pdbqt,
        display_pdb=protonated_pdb,
        cache_key=cache_key,
    )
