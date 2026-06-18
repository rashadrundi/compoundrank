import csv
import tempfile
import unittest
from pathlib import Path

from compoundrank.run_report import write_run_report


class RunReportStage4ATests(unittest.TestCase):
    def test_run_report_includes_stage4a_ligand_retrieval(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            stage4a = output_dir / "stage4a_compound_retrieval"
            ligands = stage4a / "retrieved_ligands"
            ligands.mkdir(parents=True)

            candidate_csv = stage4a / "candidate_ligands.csv"
            with candidate_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "compound_name",
                        "retrieval_rank",
                        "design_status",
                        "evidence_level",
                        "retrieval_rule_id",
                        "target_family_basis",
                        "special_domain_label",
                        "special_domain_accession",
                        "pubchem_cid",
                        "structure_fetch_status",
                        "local_sdf_path",
                        "retrieval_reason",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "compound_name": "darunavir",
                        "retrieval_rank": "1",
                        "design_status": "known_inhibitor",
                        "evidence_level": "strong",
                        "retrieval_rule_id": "retroviral_aspartyl_protease",
                        "target_family_basis": "retroviral aspartyl protease",
                        "special_domain_label": "Retroviral aspartyl protease domain",
                        "special_domain_accession": "pfam00077",
                        "pubchem_cid": "213039",
                        "structure_fetch_status": "fetched",
                        "local_sdf_path": str(ligands / "darunavir.sdf"),
                        "retrieval_reason": "Selected by rule because target evidence matched pfam00077.",
                    }
                )

            (stage4a / "docking_manifest.csv").write_text(
                "name,source_type,value\n"
                f"darunavir,file,{ligands / 'darunavir.sdf'}\n",
                encoding="utf-8",
            )
            (stage4a / "ligand_search_report.md").write_text(
                "# Ligand Search Report\n",
                encoding="utf-8",
            )

            report_path = write_run_report(output_dir=output_dir)
            text = report_path.read_text(encoding="utf-8")

            self.assertIn("## Stage 4A Ligand Retrieval", text)
            self.assertIn("darunavir", text)
            self.assertIn("retroviral_aspartyl_protease", text)
            self.assertIn("213039", text)
            self.assertIn("pfam00077", text)


if __name__ == "__main__":
    unittest.main()
