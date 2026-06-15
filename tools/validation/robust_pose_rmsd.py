from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS


def _mapped_rmsd(
    docked: Chem.Mol,
    reference: Chem.Mol,
    docked_match: tuple[int, ...],
    reference_match: tuple[int, ...],
) -> float:
    docked_conf = docked.GetConformer()
    reference_conf = reference.GetConformer()

    squared_distances = []

    for docked_index, reference_index in zip(
        docked_match,
        reference_match,
    ):
        docked_point = docked_conf.GetAtomPosition(docked_index)
        reference_point = reference_conf.GetAtomPosition(
            reference_index
        )

        squared_distances.append(
            (docked_point.x - reference_point.x) ** 2
            + (docked_point.y - reference_point.y) ** 2
            + (docked_point.z - reference_point.z) ** 2
        )

    return float(np.sqrt(np.mean(squared_distances)))


def calculate_symmetry_aware_rmsd(
    docked_molecule: Chem.Mol,
    reference_molecule: Chem.Mol,
) -> tuple[float, int]:
    docked = Chem.RemoveHs(docked_molecule)
    reference = Chem.RemoveHs(reference_molecule)

    docked_count = docked.GetNumAtoms()
    reference_count = reference.GetNumAtoms()

    if docked_count != reference_count:
        raise RuntimeError(
            "Docked and reference molecules have different "
            f"heavy-atom counts: {docked_count} versus "
            f"{reference_count}"
        )

    # First try the strict graph mapping used previously.
    strict_matches = reference.GetSubstructMatches(
        docked,
        uniquify=False,
        useChirality=True,
        maxMatches=100000,
    )

    if not strict_matches:
        strict_matches = reference.GetSubstructMatches(
            docked,
            uniquify=False,
            useChirality=False,
            maxMatches=100000,
        )

    if strict_matches:
        identity_match = tuple(range(docked_count))

        rmsds = [
            _mapped_rmsd(
                docked,
                reference,
                identity_match,
                reference_match,
            )
            for reference_match in strict_matches
        ]

        return min(rmsds), len(strict_matches)

    # GNINA refinement can produce SDF records whose bond orders or
    # formal charges differ from the reference representation. Use a
    # full-heavy-atom MCS that preserves elements and connectivity but
    # ignores bond-order differences.
    mcs = rdFMCS.FindMCS(
        [docked, reference],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        matchValences=False,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        timeout=30,
    )

    if mcs.canceled:
        raise RuntimeError("MCS atom mapping timed out")

    if mcs.numAtoms != docked_count:
        raise RuntimeError(
            "Could not establish a complete heavy-atom mapping. "
            f"MCS covered {mcs.numAtoms} of {docked_count} atoms."
        )

    query = Chem.MolFromSmarts(mcs.smartsString)

    if query is None:
        raise RuntimeError("Could not construct the MCS query")

    docked_matches = docked.GetSubstructMatches(
        query,
        uniquify=False,
        useChirality=False,
        maxMatches=100000,
    )

    reference_matches = reference.GetSubstructMatches(
        query,
        uniquify=False,
        useChirality=False,
        maxMatches=100000,
    )

    if not docked_matches or not reference_matches:
        raise RuntimeError(
            "MCS was found, but atom matches could not be enumerated"
        )

    best_rmsd = float("inf")
    mappings_tested = 0

    for docked_match in docked_matches:
        for reference_match in reference_matches:
            rmsd = _mapped_rmsd(
                docked,
                reference,
                docked_match,
                reference_match,
            )

            best_rmsd = min(best_rmsd, rmsd)
            mappings_tested += 1

    return best_rmsd, mappings_tested
