from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from compoundrank.clustering import cluster_pose_hypotheses
from compoundrank.models import PoseRecord
from compoundrank.uncertainty import assess_uncertainty


def molecule_at_offset(offset: float) -> Chem.Mol:
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(molecule, randomSeed=2026)
    molecule = Chem.RemoveHs(molecule)
    conformer = molecule.GetConformer()
    for index in range(molecule.GetNumAtoms()):
        point = conformer.GetAtomPosition(index)
        point.x += offset
        conformer.SetAtomPosition(index, point)
    return molecule


class ClusteringTests(unittest.TestCase):
    def test_clusters_use_gnina_order_and_geometry(self) -> None:
        records = [
            PoseRecord("x", 1, 1, molecule_at_offset(0.0), 0.9, 8.0, -7.0, Path("a")),
            PoseRecord("x", 2, 1, molecule_at_offset(0.2), 0.8, 8.1, -7.1, Path("b")),
            PoseRecord("x", 3, 1, molecule_at_offset(5.0), 0.7, 7.5, -6.5, Path("c")),
        ]
        clusters = cluster_pose_hypotheses(records, rmsd_threshold=1.0)
        self.assertEqual(len(clusters), 2)
        self.assertEqual(clusters[0].representative.cnn_score, 0.9)
        self.assertEqual(clusters[0].member_count, 2)
        self.assertEqual(clusters[1].member_count, 1)

    def test_uncertainty_is_descriptive(self) -> None:
        records = [
            PoseRecord("x", 1, 1, molecule_at_offset(0.0), 0.90, 8.0, -7.0, Path("a")),
            PoseRecord("x", 2, 1, molecule_at_offset(0.2), 0.88, 8.0, -7.0, Path("b")),
            PoseRecord("x", 3, 1, molecule_at_offset(5.0), 0.86, 8.0, -7.0, Path("c")),
        ]
        clusters = cluster_pose_hypotheses(records, rmsd_threshold=1.0)
        level, reasons = assess_uncertainty(clusters, seed_count=3)
        self.assertIn(level, {"low", "moderate", "high"})
        self.assertTrue(reasons)


if __name__ == "__main__":
    unittest.main()
