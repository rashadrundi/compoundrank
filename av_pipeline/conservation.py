import pandas as pd
from .utils import safe_read_csv, confidence_from_score

def parse_consurf(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: ConSurf ###
    Expected columns:
    chain,residue_number,residue_name,consurf_grade,conservation_score,exposure
    """
    df = safe_read_csv(path)
    for c in ["residue_number", "consurf_grade", "conservation_score"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def summarize_conservation(df):
    if df is None or df.empty:
        return {"conservation_score": 0, "conservation_confidence": "FAIL", "notes": "No ConSurf/custom conservation scores parsed."}
    score, notes = 0, []
    if "consurf_grade" in df.columns:
        pct = 100 * (df["consurf_grade"] >= 7).sum() / len(df)
        score += min(60, pct)
        notes.append(f"{pct:.1f}% residues grade >=7.")
    elif "conservation_score" in df.columns:
        vals = pd.to_numeric(df["conservation_score"], errors="coerce").dropna()
        if len(vals):
            mean = vals.mean()
            score += min(60, mean * 60 if mean <= 1 else mean)
    if "exposure" in df.columns and df["exposure"].astype(str).str.lower().str.contains("exposed|surface").any():
        score += 20
        notes.append("Surface/exposure labels present.")
    score += 20 if len(df) >= 50 else 10 if len(df) >= 10 else 0
    score = round(min(100, score), 2)
    return {"conservation_score": score, "conservation_confidence": confidence_from_score(score), "notes": " ".join(notes)}

def conserved_residue_set(df):
    if df is None or df.empty or "chain" not in df.columns or "residue_number" not in df.columns:
        return set()
    if "consurf_grade" in df.columns:
        sub = df[df["consurf_grade"] >= 7]
    elif "conservation_score" in df.columns:
        vals = pd.to_numeric(df["conservation_score"], errors="coerce")
        sub = df[vals >= vals.quantile(0.75)]
    else:
        return set()
    return set((str(r["chain"]), int(r["residue_number"])) for _, r in sub.dropna(subset=["residue_number"]).iterrows())
