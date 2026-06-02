from pathlib import Path
from Bio import AlignIO
import pandas as pd
import math
from .utils import confidence_from_score

def conservation_from_alignment(path, output_csv=None):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    try:
        aln = AlignIO.read(str(p), "fasta")
    except Exception:
        return pd.DataFrame()
    rows = []
    n = len(aln)
    for i in range(aln.get_alignment_length()):
        col = [str(rec.seq[i]) for rec in aln]
        non_gap = [x for x in col if x != "-"]
        if non_gap:
            counts = pd.Series(non_gap).value_counts()
            consensus = counts.index[0]
            frac = counts.iloc[0] / len(non_gap)
            probs = counts / counts.sum()
            entropy = -sum(float(p) * math.log(float(p), 2) for p in probs)
        else:
            consensus, frac, entropy = "-", 0, None
        rows.append({
            "alignment_position": i + 1,
            "consensus": consensus,
            "non_gap_count": len(non_gap),
            "num_sequences": n,
            "gap_fraction": 1 - (len(non_gap) / n if n else 0),
            "conservation_fraction": frac,
            "entropy": entropy,
        })
    df = pd.DataFrame(rows)
    if output_csv:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_csv, index=False)
    return df

def compare_msa(mafft_df, muscle_df):
    if mafft_df is None or muscle_df is None or mafft_df.empty or muscle_df.empty:
        return {"msa_comparison_score": 0, "msa_confidence": "FAIL", "notes": "MAFFT and/or MUSCLE conservation table missing."}
    n = min(len(mafft_df), len(muscle_df))
    a = mafft_df.head(n)["conservation_fraction"].astype(float)
    b = muscle_df.head(n)["conservation_fraction"].astype(float)
    agreement = max(0, 1 - (a - b).abs().mean())
    score = round(agreement * 100, 2)
    return {"msa_comparison_score": score, "msa_confidence": confidence_from_score(score, 80, 60, 35), "notes": f"MAFFT/MUSCLE conservation agreement: {score}%."}
