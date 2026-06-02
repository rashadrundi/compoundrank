import pandas as pd
from .utils import safe_read_csv, confidence_from_score
from pathlib import Path

def parse_candidates(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: ChEMBL ###
    ### EXTERNAL TOOL REQUIRED HERE: BindingDB ###
    ### EXTERNAL TOOL REQUIRED HERE: PubChem Structure Search ###
    ### EXTERNAL TOOL REQUIRED HERE: ZINC-22 ###
    Expected columns:
    ligand_id,name,smiles,source,target_class,evidence_type,evidence_notes,known_activity_value,known_activity_units
    """
    return safe_read_csv(path)

def rdkit_score(df):
    if df is None or df.empty:
        return pd.DataFrame()
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, AllChem
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        pains = FilterCatalog(params)
    except Exception as e:
        out = df.copy()
        out["ligand_score"] = 0
        out["ligand_confidence"] = "FAIL"
        out["rdkit_status"] = f"RDKit unavailable: {e}"
        return out
    rows, seen = [], set()
    for _, r in df.iterrows():
        row = r.to_dict()
        smi = str(row.get("smiles", "")).strip()
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            row.update({"rdkit_status": "invalid_smiles", "ligand_score": 0, "ligand_confidence": "FAIL"})
            rows.append(row); continue
        can = Chem.MolToSmiles(mol, canonical=True)
        duplicate = can in seen
        seen.add(can)
        mw = Descriptors.MolWt(mol); logp = Crippen.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol); hba = rdMolDescriptors.CalcNumHBA(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol); rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        alerts = [m.GetDescription() for m in pains.GetMatches(mol)]
        score, notes = 0, []
        if str(row.get("target_class","")).strip():
            score += 20; notes.append("Target-class rationale present.")
        evidence = " ".join([str(row.get("source","")), str(row.get("evidence_type","")), str(row.get("evidence_notes",""))]).lower()
        if any(k in evidence for k in ["chembl","bindingdb","bioactivity","ic50","ki","kd","ec50","inhibitor"]):
            score += 25; notes.append("Bioactivity/binding evidence language present.")
        elif any(k in evidence for k in ["pubchem","zinc","similar","analog"]):
            score += 12; notes.append("Similarity/library evidence present.")
        if not duplicate: score += 5
        lipinski = mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10
        veber = rotb <= 10 and tpsa <= 140
        if lipinski: score += 15
        if veber: score += 10
        if not alerts: score += 15
        if 150 <= mw <= 700 and rotb <= 15: score += 10
        score = min(100, score)
        row.update({"canonical_smiles": can, "duplicate": duplicate, "mw": round(mw,3), "logp": round(logp,3), "hbd": hbd, "hba": hba, "tpsa": round(tpsa,3), "rotatable_bonds": rotb, "lipinski_like": lipinski, "veber_like": veber, "pains_alerts": "; ".join(alerts), "ligand_score": score, "ligand_confidence": confidence_from_score(score), "rdkit_status": "ok", "ligand_notes": " ".join(notes)})
        rows.append(row)
    return pd.DataFrame(rows).sort_values("ligand_score", ascending=False)

def write_clean_sdf(df, output_sdf, max_ligands=50):
    """
    RDKit starter 3D SDF generation. Top ligands still need manual protonation/tautomer/stereochemistry checks.

    ### EXTERNAL TOOL REQUIRED HERE: Open Babel ###
    ### EXTERNAL TOOL REQUIRED HERE: Meeko ligand prep ###
    """
    if df is None or df.empty:
        return {"written": 0, "notes": "No ligands."}
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as e:
        return {"written": 0, "notes": f"RDKit unavailable: {e}"}
    output_sdf = Path(output_sdf); output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf)); written = 0
    for _, r in df.head(max_ligands).iterrows():
        mol = Chem.MolFromSmiles(str(r.get("canonical_smiles", r.get("smiles",""))))
        if mol is None: continue
        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            try: AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            except Exception: pass
        mol.SetProp("_Name", str(r.get("name", r.get("ligand_id", f"ligand_{written+1}"))))
        writer.write(mol); written += 1
    writer.close()
    return {"written": written, "output_sdf": str(output_sdf)}
