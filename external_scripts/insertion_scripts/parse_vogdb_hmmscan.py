import argparse
import csv
from pathlib import Path


def parse_tblout(path, max_hits=10):
    rows = []

    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return rows

    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip() or line.startswith("#"):
            continue

        parts = line.split(maxsplit=18)

        if len(parts) < 18:
            continue

        target_name = parts[0]
        target_accession = parts[1]
        query_name = parts[2]
        query_accession = parts[3]
        full_evalue = parts[4]
        full_score = parts[5]
        full_bias = parts[6]
        best_domain_evalue = parts[7]
        best_domain_score = parts[8]
        best_domain_bias = parts[9]
        exp = parts[10]
        reg = parts[11]
        clu = parts[12]
        ov = parts[13]
        env = parts[14]
        dom = parts[15]
        rep = parts[16]
        inc = parts[17]
        description = parts[18] if len(parts) > 18 else ""

        rows.append({
            "query": query_name,
            "vog_id": target_name,
            "accession": target_accession,
            "description": description,
            "evalue": full_evalue,
            "score": full_score,
            "bias": full_bias,
            "best_domain_evalue": best_domain_evalue,
            "best_domain_score": best_domain_score,
            "notes": f"hmmscan VOGDB hit; exp={exp}; domains_reported={dom}; included={inc}",
        })

    def evalue_float(row):
        try:
            return float(row["evalue"])
        except Exception:
            return 999.0

    rows.sort(key=evalue_float)
    return rows[:max_hits]


def write_csv(rows, out_path):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "query",
        "vog_id",
        "accession",
        "description",
        "evalue",
        "score",
        "bias",
        "best_domain_evalue",
        "best_domain_score",
        "notes",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tblout", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-hits", type=int, default=10)
    args = parser.parse_args()

    rows = parse_tblout(args.tblout, max_hits=args.max_hits)
    write_csv(rows, args.out)

    print(f"[VOGDB] Parsed {len(rows)} hit(s) into:")
    print(args.out)


if __name__ == "__main__":
    main()