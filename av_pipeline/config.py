from pathlib import Path
import yaml

def load_config(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not cfg or "run_dir" not in cfg:
        raise ValueError("Config must include run_dir")
    if "input" not in cfg or "protein_fasta" not in cfg["input"]:
        raise ValueError("Config must include input.protein_fasta")
    return cfg
