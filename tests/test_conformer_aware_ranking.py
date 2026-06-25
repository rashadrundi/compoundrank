from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from compoundrank.clustering import (
    cluster_pose_hypotheses,
)
from compoundrank.conformer_context import (
    receptor_display_pdb_for_pose,
)
from compoundrank.export import (
    write_complex_pdb,
)
from compoundrank.models import (
    PoseCluster,
    PoseRecord,
)


def pose(
    *,
    conformer_id: str,
    display_pdb: Path | None = None,
    seed: int = 1,
) -> PoseRecord:
    return PoseRecord(
        ligand_name="ligand",
        seed=seed,
        pose_number=1,
        molecule=object(),
        cnn_score=0.9,
        cnn_affinity=7.0,
        minimized_affinity=-8.0,
        source_sdf=Path(
            "/tmp/poses.sdf"
        ),
        receptor_conformer_id=(
            conformer_id
        ),
        receptor_display_pdb=(
            display_pdb
        ),
    )


class ConformerAwareRankingTests(
    unittest.TestCase
):
    @patch(
        "compoundrank.clustering."
        "direct_heavy_rmsd",
        return_value=0.0,
    )
    def test_different_conformers_do_not_merge(
        self,
        direct_rmsd,
    ) -> None:
        clusters = (
            cluster_pose_hypotheses(
                [
                    pose(
                        conformer_id=(
                            "snapshot_0001"
                        ),
                        seed=1,
                    ),
                    pose(
                        conformer_id=(
                            "snapshot_0002"
                        ),
                        seed=2,
                    ),
                ]
            )
        )

        self.assertEqual(
            len(clusters),
            2,
        )

        self.assertEqual(
            {
                cluster
                .representative
                .receptor_conformer_id
                for cluster in clusters
            },
            {
                "snapshot_0001",
                "snapshot_0002",
            },
        )

        direct_rmsd.assert_not_called()

    @patch(
        "compoundrank.clustering."
        "direct_heavy_rmsd",
        return_value=0.0,
    )
    def test_same_conformer_can_merge(
        self,
        direct_rmsd,
    ) -> None:
        clusters = (
            cluster_pose_hypotheses(
                [
                    pose(
                        conformer_id=(
                            "snapshot_0001"
                        ),
                        seed=1,
                    ),
                    pose(
                        conformer_id=(
                            "snapshot_0001"
                        ),
                        seed=2,
                    ),
                ]
            )
        )

        self.assertEqual(
            len(clusters),
            1,
        )

        self.assertEqual(
            clusters[0].member_count,
            2,
        )

        direct_rmsd.assert_called_once()

    def test_pose_receptor_path_or_fallback(
        self,
    ) -> None:
        conformer_path = Path(
            "/tmp/snapshot_display.pdb"
        )

        fallback = Path(
            "/tmp/submitted_display.pdb"
        )

        self.assertEqual(
            receptor_display_pdb_for_pose(
                pose(
                    conformer_id=(
                        "snapshot_0001"
                    ),
                    display_pdb=(
                        conformer_path
                    ),
                ),
                fallback,
            ),
            conformer_path,
        )

        self.assertEqual(
            receptor_display_pdb_for_pose(
                pose(
                    conformer_id=(
                        "submitted_receptor"
                    ),
                    display_pdb=None,
                ),
                fallback,
            ),
            fallback,
        )

    @patch(
        "compoundrank.export."
        "_ligand_pdb_records",
        return_value=[
            (
                "HETATM    1  C1  LIG Z   1"
                "       0.000   0.000   0.000"
            )
        ],
    )
    def test_export_records_receptor_conformer(
        self,
        ligand_records,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            receptor_path = (
                root / "receptor.pdb"
            )

            output_path = (
                root / "complex.pdb"
            )

            receptor_path.write_text(
                (
                    "ATOM      1  CA  ALA A   1"
                    "       0.000   0.000   0.000\n"
                    "END\n"
                ),
                encoding="utf-8",
            )

            record = pose(
                conformer_id=(
                    "snapshot_0001"
                ),
                display_pdb=(
                    receptor_path
                ),
            )

            cluster = PoseCluster(
                cluster_id=1,
                representative=record,
                members=[record],
                seeds={1},
                member_count=1,
                valid_member_count=1,
            )

            ligand_result = (
                SimpleNamespace(
                    ligand=SimpleNamespace(
                        name="ligand"
                    ),
                    uncertainty="low",
                    uncertainty_reasons=[],
                )
            )

            interaction_evidence = (
                SimpleNamespace(
                    closest_residue_distance=None,
                    contact_residues=[],
                    polar_contact_candidates=[],
                    hydrophobic_contact_residues=[],
                )
            )

            write_complex_pdb(
                output_path,
                receptor_path,
                ligand_result,
                cluster,
                interaction_evidence,
                compound_priority_rank=1,
                hypothesis_rank=1,
                hypothesis_count=1,
            )

            output = (
                output_path.read_text(
                    encoding="utf-8"
                )
            )

            self.assertIn(
                (
                    "RECEPTOR CONFORMER "
                    "snapshot_0001"
                ),
                output,
            )


if __name__ == "__main__":
    unittest.main()
