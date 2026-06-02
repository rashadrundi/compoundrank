import pandas as pd
from .utils import safe_read_csv, confidence_from_score

def parse_ring_edges(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: RING / Cytoscape / RINalyzer ###
    Expected columns:
    res1_chain,res1_number,res1_name,res2_chain,res2_number,res2_name,interaction_type,weight
    """
    df = safe_read_csv(path)
    for c in ["res1_number", "res2_number", "weight"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def centrality_from_edges(edges):
    if edges is None or edges.empty:
        return pd.DataFrame()
    try:
        import networkx as nx
        g = nx.Graph()
        for _, e in edges.iterrows():
            try:
                a = (str(e["res1_chain"]), int(e["res1_number"]))
                b = (str(e["res2_chain"]), int(e["res2_number"]))
                g.add_edge(a, b, weight=float(e.get("weight", 1) or 1))
            except Exception:
                pass
        deg = dict(g.degree())
        bet = nx.betweenness_centrality(g) if g.number_of_nodes() <= 2000 else {n:0 for n in g.nodes}
        return pd.DataFrame([{"chain": n[0], "residue_number": n[1], "degree": deg.get(n,0), "betweenness": bet.get(n,0)} for n in g.nodes])
    except Exception:
        rows = []
        for _, e in edges.iterrows():
            rows.append((str(e["res1_chain"]), int(e["res1_number"])))
            rows.append((str(e["res2_chain"]), int(e["res2_number"])))
        s = pd.Series(rows).value_counts()
        return pd.DataFrame([{"chain": k[0], "residue_number": k[1], "degree": int(v)} for k, v in s.items()])

def summarize_psn(centrality):
    if centrality is None or centrality.empty:
        return {"psn_score": 0, "psn_confidence": "FAIL", "notes": "No PSN/RING edges parsed."}
    score = 40 + (30 if len(centrality) >= 100 else 20 if len(centrality) >= 30 else 10)
    if "degree" in centrality.columns and centrality["degree"].max() >= 5:
        score += 20
    if "betweenness" in centrality.columns and centrality["betweenness"].max() > 0:
        score += 10
    score = min(100, score)
    return {"psn_score": score, "psn_confidence": confidence_from_score(score), "notes": f"PSN parsed for {len(centrality)} residues."}

def central_residue_set(centrality, q=0.85):
    if centrality is None or centrality.empty or "degree" not in centrality.columns:
        return set()
    cutoff = centrality["degree"].quantile(q)
    sub = centrality[centrality["degree"] >= cutoff]
    return set((str(r["chain"]), int(r["residue_number"])) for _, r in sub.iterrows())
