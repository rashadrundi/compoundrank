from pathlib import Path
import pandas as pd
from .setup_run import create_run_dirs, write_external_markers
from .fasta import summarize_fasta
from .external_tools import maybe_run_external_tools
from .annotation import parse_interpro_tsv, parse_cdd_csv, parse_pfam_tblout, normalize_hits, score_annotation
from .homologs import parse_vogdb_hits, summarize_homolog_fasta, score_homologs
from .msa import conservation_from_alignment, compare_msa
from .structure_quality import parse_plddt_from_pdb, summarize_plddt, parse_pae_summary, parse_molprobity, parse_qmean
from .conservation import parse_consurf, summarize_conservation, conserved_residue_set
from .psn import parse_ring_edges, centrality_from_edges, summarize_psn, central_residue_set
from .pockets import parse_fpocket_dir, parse_p2rank, parse_dogsite, prioritize_pockets
from .receptor_prep import load_status, score_status
from .ligands import parse_candidates, rdkit_score, write_clean_sdf
from .docking import parse_gnina, score_docking
from .scoring import make_candidates
from .report import write_outputs

def p(run_dir, rel):
    return str(Path(run_dir) / rel)

def run_pipeline(config):
    run_dir = Path(config["run_dir"])
    create_run_dirs(run_dir)
    write_external_markers(run_dir)

    print("[STEP 1] FASTA")
    fasta_summary = summarize_fasta(config["input"]["protein_fasta"])

    maybe_run_external_tools(config)

    paths = config.get("paths", {})
    print("[STEP 2-5] Annotation")
    interpro = parse_interpro_tsv(paths.get("interpro_tsv", p(run_dir, "annotation/interpro/interpro.tsv")))
    cdd = parse_cdd_csv(paths.get("cdd_csv", p(run_dir, "annotation/cdd/cdd_results.csv")))
    pfam = parse_pfam_tblout(paths.get("pfam_tbl", p(run_dir, "annotation/pfam/pfam.tbl")))
    annotation_hits = normalize_hits(interpro, cdd, pfam)
    annotation_summary = score_annotation(annotation_hits)

    print("[STEP 6-7] VOGDB/homologs")
    vogdb = parse_vogdb_hits(paths.get("vogdb_hits", p(run_dir, "homologs/vogdb_hits.csv")))
    hom_summary_raw = summarize_homolog_fasta(paths.get("homolog_fasta", p(run_dir, "homologs/homologs.fasta")))
    homolog_summary = score_homologs(vogdb, hom_summary_raw)

    print("[STEP 8] MSA")
    mafft_cons = conservation_from_alignment(paths.get("mafft_alignment", p(run_dir, "msa/mafft_alignment.fasta")), p(run_dir, "msa/mafft_conservation.csv"))
    muscle_cons = conservation_from_alignment(paths.get("muscle_alignment", p(run_dir, "msa/muscle_alignment.fasta")), p(run_dir, "msa/muscle_conservation.csv"))
    msa_summary = compare_msa(mafft_cons, muscle_cons)

    print("[STEP 9-10] Structure quality")
    plddt = parse_plddt_from_pdb(paths.get("receptor_pdb", p(run_dir, "structure/receptor.pdb")))
    structure_summary = summarize_plddt(plddt)
    pae_summary = parse_pae_summary(paths.get("pae_json", p(run_dir, "structure/pae.json")))
    molprobity_summary = parse_molprobity(paths.get("molprobity_csv", p(run_dir, "structure/molprobity.csv")))
    qmean_summary = parse_qmean(paths.get("qmean_csv", p(run_dir, "structure/qmean.csv")))

    print("[STEP 11] Conservation mapping")
    consurf = parse_consurf(paths.get("consurf_scores", p(run_dir, "conservation/consurf_scores.csv")))
    conservation_summary = summarize_conservation(consurf)
    conserved = conserved_residue_set(consurf)

    print("[STEP 12] PSN")
    ring = parse_ring_edges(paths.get("ring_edges", p(run_dir, "psn/ring_edges.csv")))
    centrality = centrality_from_edges(ring)
    psn_summary = summarize_psn(centrality)
    central = central_residue_set(centrality)

    print("[STEP 13-16] Pockets")
    fpocket = parse_fpocket_dir(paths.get("fpocket_info_dir", p(run_dir, "pockets/fpocket_out")))
    p2rank = parse_p2rank(paths.get("p2rank_csv", p(run_dir, "pockets/p2rank_predictions.csv")))
    dogsite = parse_dogsite(paths.get("dogsite_csv", p(run_dir, "pockets/dogsite_pockets.csv")))
    functional = set()
    func_csv = Path(paths.get("functional_residues_csv", p(run_dir, "conservation/functional_residues.csv")))
    if func_csv.exists():
        tmp = pd.read_csv(func_csv)
        if {"chain","residue_number"}.issubset(tmp.columns):
            functional = set((str(r["chain"]), int(r["residue_number"])) for _, r in tmp.iterrows())
    pocket_priority = prioritize_pockets(fpocket, p2rank, dogsite, conserved, central, functional)

    print("[STEP 17-19] Receptor prep")
    receptor_status = load_status(paths.get("receptor_prep_status", p(run_dir, "receptor/receptor_prep_status.yaml")))
    receptor_summary = score_status(receptor_status)

    print("[STEP 20-24] Ligands")
    ligand_raw = parse_candidates(paths.get("ligand_candidates", p(run_dir, "ligands/raw/ligand_candidates.csv")))
    ligand_priority = rdkit_score(ligand_raw)
    sdf_summary = write_clean_sdf(ligand_priority, p(run_dir, "ligands/clean/clean_ligands.sdf"))

    print("[STEP 25] Docking")
    gnina_raw = parse_gnina(paths.get("gnina_scores", p(run_dir, "docking/gnina_scores.csv")))
    docking_scored = score_docking(gnina_raw)

    print("[STEP 26-29] Candidate scoring/report")
    candidates = make_candidates(annotation_summary, homolog_summary, structure_summary, conservation_summary, psn_summary, receptor_summary, pocket_priority, ligand_priority, docking_scored)

    results = {
        "fasta_summary": fasta_summary,
        "annotation_hits": annotation_hits,
        "annotation_summary": annotation_summary,
        "homolog_summary": homolog_summary,
        "homolog_fasta_raw_summary": hom_summary_raw,
        "msa_summary": msa_summary,
        "mafft_conservation": mafft_cons,
        "muscle_conservation": muscle_cons,
        "plddt_table": plddt,
        "structure_summary": structure_summary,
        "pae_summary": pae_summary,
        "molprobity_summary": molprobity_summary,
        "qmean_summary": qmean_summary,
        "consurf_table": consurf,
        "conservation_summary": conservation_summary,
        "psn_centrality": centrality,
        "psn_summary": psn_summary,
        "pocket_priority_table": pocket_priority,
        "receptor_summary": receptor_summary,
        "ligand_priority_table": ligand_priority,
        "ligand_sdf_summary": sdf_summary,
        "docking_table": docking_scored,
        "candidate_table": candidates,
    }
    outputs = write_outputs(run_dir, results)
    print("[DONE]")
    for k, v in outputs.items():
        print(f"{k}: {v}")
    return results
