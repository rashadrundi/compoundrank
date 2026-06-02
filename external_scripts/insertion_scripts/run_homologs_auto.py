import argparse
import csv
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "antiviral-stage1-homolog-fetcher"


def read_fasta(path):
    records = []
    current_id = None
    seq_parts = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id:
                    records.append((current_id, "".join(seq_parts)))
                current_id = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    if current_id:
        records.append((current_id, "".join(seq_parts)))

    return records


def write_fasta(records, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")


def clean_header(header):
    header = header.replace("\n", " ").replace("\r", " ")
    header = re.sub(r"\s+", " ", header).strip()
    return header


def parse_fasta_text(text):
    records = []
    current_header = None
    seq_parts = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_header:
                records.append((clean_header(current_header), "".join(seq_parts)))
            current_header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)

    if current_header:
        records.append((clean_header(current_header), "".join(seq_parts)))

    return records


def infer_search_query(run_dir):
    annotation_file = Path(run_dir) / "reports" / "annotation_hits.csv"
    cdd_file = Path(run_dir) / "annotation" / "cdd" / "cdd_results.csv"

    text = ""

    for path in [annotation_file, cdd_file]:
        if path.exists():
            text += path.read_text(encoding="utf-8", errors="ignore").lower() + "\n"

    if (
        "influenza_na" in text
        or "neuraminidase" in text
        or "sialidase" in text
        or "cd15483" in text
    ):
        return '("influenza A virus"[Organism] OR "influenza B virus"[Organism]) AND neuraminidase[Title] AND srcdb_refseq[PROP]'

    # Fallback: broad viral protein query. This should rarely be used.
    return 'txid10239[Organism] AND viral protein[Title] AND srcdb_refseq[PROP]'


def eutils_get(endpoint, params, email=None, api_key=None):
    params = dict(params)
    params["tool"] = TOOL_NAME

    if email:
        params["email"] = email

    if api_key:
        params["api_key"] = api_key

    url = f"{EUTILS_BASE}/{endpoint}?" + urllib.parse.urlencode(params)

    req = urllib.request.Request(
        url,
        headers={"User-Agent": f"{TOOL_NAME}/0.1 educational defensive-use"},
    )

    with urllib.request.urlopen(req, timeout=120) as response:
        return response.read().decode("utf-8", errors="replace")


def esearch_protein(term, retmax, email=None, api_key=None):
    params = {
        "db": "protein",
        "term": term,
        "retmax": str(retmax),
        "retmode": "json",
        "sort": "relevance",
    }

    text = eutils_get("esearch.fcgi", params, email=email, api_key=api_key)
    data = json.loads(text)
    return data.get("esearchresult", {}).get("idlist", []), data


def efetch_fasta(ids, email=None, api_key=None):
    if not ids:
        return ""

    params = {
        "db": "protein",
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
    }

    return eutils_get("efetch.fcgi", params, email=email, api_key=api_key)


def aa_clean(seq):
    seq = seq.upper().replace("*", "")
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWYXBZUOJ\-]", "", seq)
    return seq


def dedupe_records(records):
    seen_sequences = set()
    clean_records = []

    for header, seq in records:
        seq = aa_clean(seq)

        if not seq:
            continue

        if seq in seen_sequences:
            continue

        seen_sequences.add(seq)
        clean_records.append((header, seq))

    return clean_records


def write_search_metadata(run_dir, query, ids, esearch_json):
    out_dir = Path(run_dir) / "homologs"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "homolog_search_query.txt", "w", encoding="utf-8") as f:
        f.write(query + "\n")

    with open(out_dir / "homolog_ncbi_ids.txt", "w", encoding="utf-8") as f:
        for protein_id in ids:
            f.write(str(protein_id) + "\n")

    with open(out_dir / "homolog_esearch_response.json", "w", encoding="utf-8") as f:
        json.dump(esearch_json, f, indent=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--max-records", type=int, default=25)
    parser.add_argument("--email", default="")
    parser.add_argument("--api-key", default="")
    parser.add_argument("--query", default="")
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    input_fasta = run_dir / "input" / "protein.fasta"
    output_fasta = run_dir / "homologs" / "homologs.fasta"

    if not input_fasta.exists():
        raise FileNotFoundError(f"Missing input FASTA: {input_fasta}")

    query_records = read_fasta(input_fasta)

    if not query_records:
        raise RuntimeError(f"No FASTA records found in: {input_fasta}")

    query_id, query_seq = query_records[0]

    search_query = args.query.strip() or infer_search_query(run_dir)

    print("[Homologs] Search query:")
    print(search_query)

    ids, esearch_json = esearch_protein(
        search_query,
        retmax=args.max_records,
        email=args.email,
        api_key=args.api_key,
    )

    print(f"[Homologs] Found {len(ids)} NCBI Protein IDs.")

    write_search_metadata(run_dir, search_query, ids, esearch_json)

    time.sleep(1)

    fasta_text = efetch_fasta(
        ids,
        email=args.email,
        api_key=args.api_key,
    )

    fetched_records = parse_fasta_text(fasta_text)

    print(f"[Homologs] Downloaded {len(fetched_records)} FASTA records.")

    all_records = [(f"QUERY__{query_id}", aa_clean(query_seq))] + fetched_records
    all_records = dedupe_records(all_records)

    write_fasta(all_records, str(output_fasta))

    print(f"[Homologs] Wrote {len(all_records)} total records to:")
    print(output_fasta)


if __name__ == "__main__":
    main()