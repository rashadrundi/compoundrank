from __future__ import annotations

import json
import shutil
import subprocess
import time
import traceback
from pathlib import Path
from typing import Any

from .target_evidence import build_target_evidence, write_target_evidence_outputs

DEFAULT_API_URL = "http://161.35.0.191:8000/analyze/fasta"


def post_fasta(
    api_url: str,
    fasta_path: Path,
    timeout_seconds: int = 7200,
) -> dict[str, Any]:
    fasta_path = Path(fasta_path).resolve()

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file does not exist: {fasta_path}")

    curl_path = shutil.which("curl")

    if curl_path is None:
        raise RuntimeError("curl was not found on PATH")

    command = [
        curl_path,
        "--silent",
        "--show-error",
        "--fail-with-body",
        "--connect-timeout",
        "20",
        "--max-time",
        str(timeout_seconds),
        "--request",
        "POST",
        api_url,
        "--header",
        "Accept: application/json",
        "--form",
        f"file=@{fasta_path};type=application/octet-stream",
    ]

    print(f"[CPU] Sending POST request to: {api_url}", flush=True)
    print(f"[CPU] FASTA: {fasta_path}", flush=True)

    started_at = time.monotonic()

    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )

    elapsed = time.monotonic() - started_at

    print(
        f"[CPU] curl finished after {elapsed:.1f} seconds "
        f"with exit code {completed.returncode}",
        flush=True,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "CPU API curl request failed.\n\n"
            f"Exit code: {completed.returncode}\n"
            f"STDERR:\n{completed.stderr}\n"
            f"RESPONSE BODY:\n{completed.stdout[:3000]}"
        )

    try:
        response_data = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            "CPU API returned a response, but it was not valid JSON.\n\n"
            f"Response:\n{completed.stdout[:3000]}\n"
            f"STDERR:\n{completed.stderr}"
        ) from error

    if not isinstance(response_data, dict):
        raise RuntimeError(
            "CPU API returned valid JSON, but the top-level value "
            "was not an object."
        )

    return response_data


def get_rows(
    api_response: dict[str, Any],
    tool_name: str,
) -> list[dict[str, Any]]:
    results = api_response.get("results", {})
    rows = results.get(tool_name)

    if isinstance(rows, list):
        return rows

    tool_object = api_response.get(tool_name, {})
    fallback_rows = tool_object.get("rows", [])

    if isinstance(fallback_rows, list):
        return fallback_rows

    return []


def parse_cpu_response(api_response: dict[str, Any]) -> dict[str, Any]:
    if isinstance(api_response.get("rows"), dict):
        rows = api_response.get("rows", {})

        normalized_rows: dict[str, list[dict[str, Any]]] = {}
        for tool_name in ("cdd", "interpro", "vogdb"):
            tool_rows = rows.get(tool_name, [])
            normalized_rows[tool_name] = tool_rows if isinstance(tool_rows, list) else []

        return {
            "job_id": api_response.get("job_id"),
            "status": api_response.get("status"),
            "files": api_response.get("files", {}),
            "result_counts": {
                "cdd": len(normalized_rows["cdd"]),
                "interpro": len(normalized_rows["interpro"]),
                "vogdb": len(normalized_rows["vogdb"]),
            },
            "rows": normalized_rows,
        }

    cdd_rows = get_rows(api_response, "cdd")
    interpro_rows = get_rows(api_response, "interpro")
    vogdb_rows = get_rows(api_response, "vogdb")

    return {
        "job_id": api_response.get("job_id"),
        "status": api_response.get("status"),
        "result_counts": {
            "cdd": len(cdd_rows),
            "interpro": len(interpro_rows),
            "vogdb": len(vogdb_rows),
        },
        "rows": {
            "cdd": cdd_rows,
            "interpro": interpro_rows,
            "vogdb": vogdb_rows,
        },
        "files": api_response.get("files", {}),
    }


def _write_json(
    path: Path,
    payload: dict[str, Any],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
            default=str,
        )
        + "\n",
        encoding="utf-8",
    )


def run_homolog_search(
    *,
    fasta_path: Path,
    output_dir: Path,
    api_url: str = DEFAULT_API_URL,
    timeout_seconds: int = 7200,
) -> dict[str, Any]:
    output_dir = Path(output_dir)
    fasta_path = Path(fasta_path)

    raw_output = output_dir / "homolog_search_raw.json"
    summary_output = output_dir / "homolog_search_summary.json"
    error_output = output_dir / "homolog_search_error.json"

    try:
        raw_response = post_fasta(
            api_url=api_url,
            fasta_path=fasta_path,
            timeout_seconds=timeout_seconds,
        )

        parsed = parse_cpu_response(raw_response)

        _write_json(
            raw_output,
            raw_response,
        )

        _write_json(
            summary_output,
            parsed,
        )

        target_evidence = build_target_evidence(
            parsed,
            source_fasta=str(fasta_path.resolve()),
        )

        target_outputs = write_target_evidence_outputs(
            target_evidence,
            output_dir,
        )

        return {
            "status": "ok",
            "api_url": api_url,
            "fasta_path": str(fasta_path.resolve()),
            "raw_output": str(raw_output),
            "summary_output": str(summary_output),
            "target_evidence": target_outputs.get("target_evidence"),
            "target_evidence_report": target_outputs.get("target_evidence_report"),
            "result_counts": parsed.get("result_counts", {}),
        }

    except Exception as error:
        payload = {
            "status": "error",
            "api_url": api_url,
            "fasta_path": str(fasta_path.resolve()),
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
        }

        _write_json(
            error_output,
            payload,
        )

        return {
            "status": "error",
            "api_url": api_url,
            "fasta_path": str(fasta_path.resolve()),
            "error_output": str(error_output),
            "error_type": type(error).__name__,
            "error": str(error),
        }
