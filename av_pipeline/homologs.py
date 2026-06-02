from pathlib import Path
import pandas as pd
from .utils import safe_read_csv


def parse_vogdb_hits(path):
    return safe_read_csv(path)


def _read_fasta_records(path):
    p = Path(path)

    if not p.exists() or p.stat().st_size == 0:
        return []

    records = []
    current_header = None
    seq_parts = []

    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(seq_parts)))
            current_header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)

    if current_header is not None:
        records.append((current_header, "".join(seq_parts)))

    return records


def summarize_homolog_fasta(path):
    records = _read_fasta_records(path)

    if not records:
        return {
            "num_sequences": 0,
            "num_homologs": 0,
            "headers": [],
            "mean_length": None,
            "notes": "No homolog FASTA records found.",
        }

    lengths = [len(seq) for _, seq in records if seq]
    headers = [header for header, _ in records]

    query_count = sum(1 for h in headers if h.upper().startswith("QUERY__"))
    num_sequences = len(records)
    num_homologs = max(0, num_sequences - query_count)

    mean_length = sum(lengths) / len(lengths) if lengths else None

    return {
        "num_sequences": num_sequences,
        "num_homologs": num_homologs,
        "headers": headers[:10],
        "mean_length": mean_length,
        "notes": f"Parsed {num_sequences} total FASTA records, including {num_homologs} non-query homolog records.",
    }


def score_homologs(vogdb, summary):
    score = 0
    notes = []

    num_sequences = 0
    num_homologs = 0

    if isinstance(summary, dict):
        num_sequences = summary.get("num_sequences", 0) or 0
        num_homologs = summary.get("num_homologs", 0) or 0

    try:
        num_sequences = int(num_sequences)
    except Exception:
        num_sequences = 0

    try:
        num_homologs = int(num_homologs)
    except Exception:
        num_homologs = 0

    if num_homologs >= 25:
        score += 50
        notes.append(f"Homolog FASTA contains {num_homologs} non-query homologs, sufficient for a strong POC MSA.")
    elif num_homologs >= 10:
        score += 40
        notes.append(f"Homolog FASTA contains {num_homologs} non-query homologs, sufficient for a useful POC MSA.")
    elif num_homologs >= 5:
        score += 30
        notes.append(f"Homolog FASTA contains {num_homologs} non-query homologs, enough for preliminary homolog context.")
    elif num_homologs >= 2:
        score += 15
        notes.append(f"Homolog FASTA contains {num_homologs} non-query homologs, but more homologs are recommended.")
    else:
        notes.append("Insufficient homolog sequences for robust MSA conservation.")

    if isinstance(summary, dict) and summary.get("mean_length"):
        notes.append(f"Mean homolog FASTA sequence length: {summary.get('mean_length'):.1f} aa.")

    if vogdb is not None and not vogdb.empty:
        score += 25
        notes.append("VOGDB/homolog-context table present.")

        text = " ".join(vogdb.astype(str).fillna("").values.flatten()).lower()

        if any(x in text for x in ["vog", "vfam", "vfold"]):
            score += 15
            notes.append("VOG/VFAM/VFOLD-style viral ortholog identifiers present.")

        if any(x in text for x in ["neuraminidase", "sialidase", "influenza"]):
            score += 15
            notes.append("VOGDB/homolog context supports influenza neuraminidase identity.")
    else:
        notes.append("No VOGDB hit table available yet; homolog confidence is based on FASTA sequence set only.")

    score = min(100, score)

    if score >= 75:
        confidence = "HIGH"
    elif score >= 50:
        confidence = "MODERATE"
    elif score >= 25:
        confidence = "LOW"
    else:
        confidence = "FAIL"

    return {
        "homolog_score": score,
        "homolog_confidence": confidence,
        "notes": " ".join(notes),
    }