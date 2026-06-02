from pathlib import Path
from .utils import ensure_dir

DIRS = [
    "input",
    "annotation/interpro",
    "annotation/cdd",
    "annotation/pfam",
    "homologs",
    "msa",
    "structure",
    "conservation",
    "psn",
    "pockets/fpocket_out",
    "receptor",
    "ligands/raw",
    "ligands/clean",
    "ligands/prepared",
    "docking",
    "reports",
    "external_markers",
    "logs",
]

MARKERS = {
    "01_interproscan.txt": "### EXTERNAL TOOL REQUIRED HERE: InterProScan ###\nPlace TSV at annotation/interpro/interpro.tsv\n",
    "02_cdd_cdsearch.txt": "### EXTERNAL TOOL REQUIRED HERE: CDD / CD-Search ###\nPlace CSV at annotation/cdd/cdd_results.csv\n",
    "03_pfam_hmmer.txt": "### EXTERNAL TOOL REQUIRED HERE: Pfam / HMMER ###\nPlace hmmscan --tblout at annotation/pfam/pfam.tbl\n",
    "04_vogdb.txt": "### EXTERNAL TOOL REQUIRED HERE: VOGDB ###\nPlace VOGDB hits at homologs/vogdb_hits.csv and homologs at homologs/homologs.fasta\n",
    "05_mafft.txt": "### EXTERNAL TOOL REQUIRED HERE: MAFFT ###\nPlace alignment at msa/mafft_alignment.fasta\n",
    "06_muscle.txt": "### EXTERNAL TOOL REQUIRED HERE: MUSCLE ###\nPlace alignment at msa/muscle_alignment.fasta\n",
    "07_alphafold_colabfold.txt": "### EXTERNAL TOOL REQUIRED HERE: AlphaFold / ColabFold ###\nPlace model at structure/receptor.pdb and optional PAE at structure/pae.json\n",
    "08_consurf.txt": "### EXTERNAL TOOL REQUIRED HERE: ConSurf ###\nPlace scores at conservation/consurf_scores.csv\n",
    "09_ring_cytoscape_rinalyzer.txt": "### EXTERNAL TOOL REQUIRED HERE: RING / Cytoscape / RINalyzer ###\nPlace edge table at psn/ring_edges.csv\n",
    "10_pockets.txt": "### EXTERNAL TOOL REQUIRED HERE: fpocket / P2Rank / PrankWeb / DoGSiteScorer ###\nPlace outputs in pockets/\n",
    "11_pdbfixer.txt": "### EXTERNAL TOOL REQUIRED HERE: PDBFixer ###\nRecord repair status in receptor/receptor_prep_status.yaml\n",
    "12_chimerax.txt": "### EXTERNAL TOOL REQUIRED HERE: ChimeraX ###\nRecord inspection status in receptor/receptor_prep_status.yaml\n",
    "13_meeko_receptor.txt": "### EXTERNAL TOOL REQUIRED HERE: Meeko receptor prep ###\nRecord receptor/box status in receptor/receptor_prep_status.yaml\n",
    "14_ligand_sources.txt": "### EXTERNAL TOOL REQUIRED HERE: ChEMBL / BindingDB / PubChem Structure Search / ZINC-22 ###\nPlace ligand candidates at ligands/raw/ligand_candidates.csv\n",
    "15_ligand_prep.txt": "### EXTERNAL TOOL REQUIRED HERE: Open Babel / Meeko ligand prep ###\nPrepare ligands after RDKit cleaning\n",
    "16_gnina.txt": "### EXTERNAL TOOL REQUIRED HERE: GNINA ###\nPlace scores at docking/gnina_scores.csv\n",
}

def create_run_dirs(run_dir):
    run_dir = Path(run_dir)
    for d in DIRS:
        ensure_dir(run_dir / d)

def write_external_markers(run_dir):
    d = ensure_dir(Path(run_dir) / "external_markers")
    for name, text in MARKERS.items():
        (d / name).write_text(text, encoding="utf-8")
