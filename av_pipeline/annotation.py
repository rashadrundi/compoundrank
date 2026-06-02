from pathlib import Path
import pandas as pd
from .utils import safe_read_csv, confidence_from_score

def parse_interpro_tsv(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    cols = [
        "protein_accession", "sequence_md5", "sequence_length", "analysis",
        "signature_accession", "signature_description", "start", "end",
        "score", "status", "date", "interpro_accession",
        "interpro_description", "go_terms", "pathways"
    ]
    df = pd.read_csv(p, sep="\t", names=cols, comment="#", dtype=str)
    for c in ["sequence_length", "start", "end"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def parse_cdd_csv(path):
    df = safe_read_csv(path)

    if df.empty:
        return df

    # Case 1:
    # Already-normalized CDD CSV from our run_cdd_auto.py script:
    # source,name,accession,start,end,evalue,score,notes
    if "name" in df.columns and "accession" in df.columns:
        out = df.copy()

        if "score" in out.columns and "bitscore" not in out.columns:
            out["bitscore"] = out["score"]

        if "accession" in out.columns and "hit_id" not in out.columns:
            out["hit_id"] = out["accession"]

        if "name" in out.columns and "short_name" not in out.columns:
            out["short_name"] = out["name"]

        if "description" not in out.columns:
            out["description"] = out.get("notes", "")

        if "superfamily" not in out.columns:
            out["superfamily"] = ""

        if "sites" not in out.columns:
            out["sites"] = ""

        if "query" not in out.columns:
            out["query"] = ""

        return out

    # Case 2:
    # Older/raw CDD-style CSV expected by the original pipeline:
    # query,hit_id,short_name,description,start,end,evalue,bitscore,superfamily,sites
    rename_map = {
        "From": "start",
        "To": "end",
        "E-Value": "evalue",
        "Bitscore": "bitscore",
        "Accession": "hit_id",
        "Short name": "short_name",
        "Superfamily": "superfamily",
    }

    df = df.rename(columns=rename_map)

    for col in [
        "query",
        "hit_id",
        "short_name",
        "description",
        "start",
        "end",
        "evalue",
        "bitscore",
        "superfamily",
        "sites",
    ]:
        if col not in df.columns:
            df[col] = ""

    return df

def parse_pfam_tblout(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: Pfam / HMMER ###
    Parses hmmscan --tblout output.
    """
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    rows = []
    for line in p.read_text(errors="ignore").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split(maxsplit=18)
        if len(parts) < 7:
            continue
        def f(x):
            try: return float(x)
            except Exception: return None
        rows.append({
            "target_name": parts[0],
            "target_accession": parts[1] if len(parts) > 1 else "",
            "query_name": parts[2] if len(parts) > 2 else "",
            "full_evalue": f(parts[4]) if len(parts) > 4 else None,
            "full_score": f(parts[5]) if len(parts) > 5 else None,
            "description": parts[18] if len(parts) > 18 else "",
        })
    return pd.DataFrame(rows)

def normalize_hits(interpro, cdd, pfam):
    rows = []
    if interpro is not None and not interpro.empty:
        for _, r in interpro.iterrows():
            rows.append({
                "source": "InterPro",
                "name": r.get("interpro_description") or r.get("signature_description"),
                "accession": r.get("interpro_accession") or r.get("signature_accession"),
                "start": r.get("start"),
                "end": r.get("end"),
                "evalue": None,
                "score": r.get("score"),
                "notes": r.get("analysis", ""),
            })
    if cdd is not None and not cdd.empty:
        for _, r in cdd.iterrows():
            rows.append({
                "source": "CDD",
                "name": r.get("short_name", r.get("description", "")),
                "accession": r.get("hit_id", ""),
                "start": r.get("start", ""),
                "end": r.get("end", ""),
                "evalue": r.get("evalue", None),
                "score": r.get("bitscore", None),
                "notes": r.get("notes") or r.get("description") or f"superfamily={r.get('superfamily','')}; sites={r.get('sites','')}",
            })
    if pfam is not None and not pfam.empty:
        for _, r in pfam.iterrows():
            rows.append({
                "source": "Pfam",
                "name": r.get("target_name", ""),
                "accession": r.get("target_accession", ""),
                "start": "",
                "end": "",
                "evalue": r.get("full_evalue", None),
                "score": r.get("full_score", None),
                "notes": r.get("description", ""),
            })
    return pd.DataFrame(rows)

def target_bucket(name):
    s = str(name).lower()
    buckets = {
    "neuraminidase": ["neuraminidase", "sialidase", "influenza_na"],
    "protease": ["protease", "peptidase", "proteinase"],
    "polymerase": ["polymerase", "rdrp", "rna-dependent", "replicase"],
    "helicase_atpase": ["helicase", "atpase", "ntpase"],
    "methyltransferase": ["methyltransferase", "mtase", "sam-dependent"],
    "nuclease": ["nuclease", "endoribonuclease", "exonuclease", "rnase", "dnase"],
    "entry_fusion": ["glycoprotein", "fusion", "spike", "envelope"],
    "structural": ["capsid", "nucleocapsid", "matrix", "tegument"],
}
    for k, words in buckets.items():
        if any(w in s for w in words):
            return k
    return "other"

def score_annotation(hits):
    if hits is None or hits.empty:
        return {"annotation_score": 0, "annotation_confidence": "FAIL", "working_annotation": "No confident annotation", "notes": "No annotation hits parsed."}

    df = hits.copy()
    df["bucket"] = df["name"].map(target_bucket)
    score = 0
    notes = []
    n_sources = df["source"].nunique()
    score += {1: 10, 2: 25}.get(n_sources, 35 if n_sources >= 3 else 0)

    useful = df[df["bucket"] != "other"]
    if not useful.empty:
        top = useful["bucket"].value_counts().index[0]
        count = useful["bucket"].value_counts().iloc[0]
        if count >= 3:
            score += 35
        elif count == 2:
            score += 25
        else:
            score += 15
        notes.append(f"Working target-class bucket: {top}")
        working = f"Predicted {top.replace('_',' ')}-related protein/domain"
    else:
        working = f"Putative annotation: {df.iloc[0].get('name','unclear')}"

    ev = pd.to_numeric(df.get("evalue", pd.Series(dtype=float)), errors="coerce")
    if len(ev.dropna()):
        if (ev < 1e-10).any():
            score += 15
            notes.append("Very strong E-value hit present.")
        elif (ev < 1e-3).any():
            score += 8
            notes.append("Moderate E-value hit present.")

    all_notes = " ".join(df.get("notes", pd.Series(dtype=str)).fillna("").astype(str)).lower()
    if any(w in all_notes for w in ["active", "catalytic", "binding", "site", "metal"]):
        score += 10
        notes.append("Functional-site language detected.")

    score = min(100, score)
    return {
        "annotation_score": score,
        "annotation_confidence": confidence_from_score(score),
        "working_annotation": working,
        "supporting_sources": sorted(df["source"].dropna().unique().tolist()),
        "notes": " ".join(notes),
    }
