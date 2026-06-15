from __future__ import annotations

import csv
import shutil
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from rdkit import Chem
from rdkit.Chem import AllChem

from .models import PreparedLigand
from .paths import content_cache_key, sanitize_name
from .subprocess_utils import resolve_executable, run_command


@dataclass(frozen=True)
class LigandRequest:
    name: str
    source_type: str
    value: str


def _download_pubchem_sdf(cid: str, destination: Path) -> str:
    urls = [
        (
            "3d",
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d",
        ),
        (
            "2d",
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=2d",
        ),
    ]
    errors: list[str] = []
    for label, url in urls:
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "compoundrank-local/0.2"},
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
            if data:
                destination.write_bytes(data)
                return label
        except (urllib.error.URLError, TimeoutError) as error:
            errors.append(f"{label}: {error}")
    raise RuntimeError(
        f"Could not retrieve PubChem CID {cid}. Attempts: {'; '.join(errors)}"
    )


def _sdf_has_3d_coordinates(path: Path) -> bool:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    molecule = next((mol for mol in supplier if mol is not None), None)
    if molecule is None or molecule.GetNumConformers() == 0:
        return False
    conformer = molecule.GetConformer()
    return any(
        abs(conformer.GetAtomPosition(index).z) > 1e-3
        for index in range(molecule.GetNumAtoms())
    )


def _write_smiles_3d(smiles: str, destination: Path, seed: int = 2026) -> None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    parameters.enforceChirality = True
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise RuntimeError("RDKit could not generate a 3D ligand conformer")
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        AllChem.MMFFOptimizeMolecule(molecule, maxIters=2000)
    elif AllChem.UFFHasAllMoleculeParams(molecule):
        AllChem.UFFOptimizeMolecule(molecule, maxIters=2000)
    writer = Chem.SDWriter(str(destination))
    writer.write(molecule)
    writer.close()


def read_manifest(path: Path) -> list[LigandRequest]:
    requests: list[LigandRequest] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"name", "source_type", "value"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError(
                "Ligand manifest requires columns: name,source_type,value"
            )
        for row in reader:
            requests.append(
                LigandRequest(
                    name=row["name"].strip(),
                    source_type=row["source_type"].strip().lower(),
                    value=row["value"].strip(),
                )
            )
    return requests


def prepare_ligand(
    request: LigandRequest,
    cache_root: Path,
    *,
    ph: float = 7.4,
    obabel_bin: str = "obabel",
    meeko_ligand_bin: str = "mk_prepare_ligand.py",
) -> PreparedLigand:
    name = sanitize_name(request.name)
    source_type = request.source_type.lower()
    if source_type not in {"file", "cid", "smiles"}:
        raise ValueError(f"Unsupported ligand source_type: {source_type}")

    if source_type == "file":
        source_file = Path(request.value).expanduser().resolve()
        if not source_file.is_file():
            raise FileNotFoundError(f"Ligand file does not exist: {source_file}")
        key_part: str | Path = source_file
    else:
        key_part = request.value

    cache_key = content_cache_key(
        str(source_type),
        key_part,
        f"ph={ph}",
        "ligand-v2",
    )
    cache_dir = cache_root / "ligands" / cache_key
    normalized_sdf = cache_dir / "ligand_normalized.sdf"
    prepared_pdbqt = cache_dir / "ligand_prepared.pdbqt"

    if normalized_sdf.is_file() and prepared_pdbqt.is_file():
        return PreparedLigand(
            name=name,
            source_description=f"{source_type}:{request.value}",
            source_sdf=normalized_sdf,
            prepared_pdbqt=prepared_pdbqt,
            cache_key=cache_key,
        )

    cache_dir.mkdir(parents=True, exist_ok=True)
    raw_sdf = cache_dir / "ligand_raw.sdf"

    if source_type == "file":
        source_file = Path(request.value).expanduser().resolve()
        if source_file.suffix.lower() in {".sdf", ".sd"}:
            shutil.copy2(source_file, raw_sdf)
        else:
            obabel = resolve_executable(obabel_bin, "Open Babel")
            run_command([obabel, str(source_file), "-O", str(raw_sdf)])
    elif source_type == "cid":
        _download_pubchem_sdf(request.value, raw_sdf)
    else:
        _write_smiles_3d(request.value, raw_sdf)

    obabel = resolve_executable(obabel_bin, "Open Babel")
    command = [
        obabel,
        str(raw_sdf),
        "-O",
        str(normalized_sdf),
        "-p",
        str(ph),
    ]
    if not _sdf_has_3d_coordinates(raw_sdf):
        command.append("--gen3d")
    run_command(command)

    meeko = resolve_executable(meeko_ligand_bin, "Meeko ligand preparation")
    run_command(
        [
            meeko,
            "-i",
            str(normalized_sdf),
            "-o",
            str(prepared_pdbqt),
            "--charge_model",
            "gasteiger",
            "--add_index_map",
        ]
    )

    if not prepared_pdbqt.is_file() or prepared_pdbqt.stat().st_size == 0:
        raise RuntimeError("Meeko did not create a ligand PDBQT")

    return PreparedLigand(
        name=name,
        source_description=f"{source_type}:{request.value}",
        source_sdf=normalized_sdf,
        prepared_pdbqt=prepared_pdbqt,
        cache_key=cache_key,
    )


def combine_requests(
    ligand_files: Iterable[str],
    ligand_cids: Iterable[str],
    ligand_smiles: Iterable[str],
    manifest_requests: Iterable[LigandRequest],
) -> list[LigandRequest]:
    requests = list(manifest_requests)
    for value in ligand_files:
        path = Path(value)
        requests.append(LigandRequest(path.stem, "file", value))
    for value in ligand_cids:
        requests.append(LigandRequest(f"pubchem_{value}", "cid", value))
    for value in ligand_smiles:
        if "=" in value:
            name, smiles = value.split("=", 1)
        else:
            name, smiles = "smiles_ligand", value
        requests.append(LigandRequest(name, "smiles", smiles))
    if not requests:
        raise ValueError(
            "Supply at least one --ligand-file, --ligand-cid, "
            "--ligand-smiles, or --ligand-manifest"
        )
    return requests
