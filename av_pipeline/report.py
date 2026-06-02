from pathlib import Path
import pandas as pd
from .utils import ensure_parent, write_json, df_records, now_iso

def md_table(df, max_rows=20):
    if df is None or df.empty:
        return "_No data available._"
    return df.head(max_rows).fillna("").to_markdown(index=False)

def write_outputs(run_dir, results):
    rd = Path(run_dir); reports = rd / "reports"; reports.mkdir(parents=True, exist_ok=True)
    for name in ["annotation_hits","mafft_conservation","muscle_conservation","plddt_table","consurf_table","psn_centrality","pocket_priority_table","ligand_priority_table","docking_table","candidate_table"]:
        df = results.get(name)
        if isinstance(df, pd.DataFrame) and not df.empty:
            df.to_csv(reports / f"{name}.csv", index=False)
    report = reports / "final_handoff_report.md"
    manifest = reports / "biophysics_handoff_manifest.json"
    write_report(report, results)
    write_manifest(manifest, results)
    return {"report": str(report), "manifest": str(manifest)}

def write_report(path, r):
    lines = []
    lines.append("# Biophysics Simulation Handoff Report\n")
    lines.append(f"Generated: {now_iso()}\n")
    lines.append("## Defensive-use statement\n")
    lines.append("This report supports therapeutic antiviral prioritization only. It is not proof of antiviral activity and must not be used for viral enhancement or gain-of-function purposes.\n")

    for title, key in [
        ("Input FASTA summary", "fasta_summary"),
        ("Annotation summary", "annotation_summary"),
        ("Homolog/VOGDB summary", "homolog_summary"),
        ("MSA comparison", "msa_summary"),
        ("Structure quality summary", "structure_summary"),
        ("PAE summary", "pae_summary"),
        ("MolProbity summary", "molprobity_summary"),
        ("QMEAN summary", "qmean_summary"),
        ("Conservation summary", "conservation_summary"),
        ("PSN summary", "psn_summary"),
        ("Receptor prep summary", "receptor_summary"),
    ]:
        lines.append(f"## {title}\n")
        val = r.get(key, {})
        if isinstance(val, dict):
            for k, v in val.items(): lines.append(f"- **{k}:** {v}")
        else:
            lines.append(str(val))
        lines.append("")

    lines.append("## Annotation hits\n" + md_table(r.get("annotation_hits")) + "\n")
    lines.append("## Pocket prioritization\n" + md_table(r.get("pocket_priority_table"), 15) + "\n")
    lines.append("## Ligand prioritization\n" + md_table(r.get("ligand_priority_table"), 15) + "\n")
    lines.append("## Docking summary\n" + md_table(r.get("docking_table"), 15) + "\n")
    lines.append("## Final candidate ranking\n" + md_table(r.get("candidate_table"), 20) + "\n")
    lines.append("## Simulation question\n")
    lines.append("For the top candidate(s), test whether the proposed ligand remains stably and plausibly bound to the prioritized conserved pocket, whether key contacts persist, and whether the result is consistent with annotation/conservation evidence.\n")
    lines.append("## Major uncertainty flags\n")
    lines.append("- Docking is hypothesis generation, not proof.\n- External-tool outputs must be reviewed for quality.\n- Ligand protonation, tautomer state, stereochemistry, and receptor prep require manual review for top candidates.\n")
    ensure_parent(path).write_text("\n".join(lines), encoding="utf-8")

def write_manifest(path, r):
    obj = {
        "generated_at": now_iso(),
        "purpose": "defensive antiviral protein-pocket-ligand prioritization for biophysics simulation handoff",
        "confidence_flags": {
            "HIGH": "Multiple evidence layers agree.",
            "MODERATE": "Useful hypothesis with missing/incomplete evidence.",
            "LOW": "Exploratory only.",
            "FAIL": "Not ready for handoff."
        },
        "summaries": {k:v for k,v in r.items() if isinstance(v, dict)},
        "top_pockets": df_records(r.get("pocket_priority_table"), 10),
        "top_ligands": df_records(r.get("ligand_priority_table"), 10),
        "top_docking": df_records(r.get("docking_table"), 10),
        "top_candidates": df_records(r.get("candidate_table"), 10),
    }
    write_json(path, obj)
