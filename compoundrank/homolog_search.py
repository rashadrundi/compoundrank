from __future__ import annotations

import json
import shutil
import subprocess
import time
import traceback
from pathlib import Path
from typing import Any

from .target_evidence import (
    build_target_evidence,
    write_target_evidence_outputs,
)

DEFAULT_API_URL = "http://161.35.0.191:8000/analyze/fasta"

TOOL_NAMES = ("cdd", "interpro", "vogdb")

SUCCESSFUL_TOOL_STATUSES = {
    "complete",
    "complete_no_hits",
}

KNOWN_TOOL_STATUSES = {
    "complete",
    "complete_no_hits",
    "partial",
    "failed",
    "skipped",
    "unknown",
}


def post_fasta(
    api_url: str,
    fasta_path: Path,
    timeout_seconds: int = 7200,
) -> dict[str, Any]:
    fasta_path = Path(fasta_path).resolve()

    if not fasta_path.exists():
        raise FileNotFoundError(
            f"FASTA file does not exist: {fasta_path}"
        )

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

    print(
        f"[CPU] Sending POST request to: {api_url}",
        flush=True,
    )
    print(
        f"[CPU] FASTA: {fasta_path}",
        flush=True,
    )

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


def _normalize_status_value(
    value: Any,
) -> str:
    if value is None:
        return "unknown"

    normalized = str(value).strip().lower()

    aliases = {
        "ok": "complete",
        "success": "complete",
        "successful": "complete",
        "completed": "complete",
        "no_hits": "complete_no_hits",
        "no_match": "complete_no_hits",
        "no_matches": "complete_no_hits",
        "empty": "complete_no_hits",
        "partial_or_failed": "partial",
        "error": "failed",
        "failure": "failed",
        "not_run": "skipped",
        "disabled": "skipped",
    }

    normalized = aliases.get(
        normalized,
        normalized,
    )

    if normalized in KNOWN_TOOL_STATUSES:
        return normalized

    return "unknown"


def _get_tool_object(
    api_response: dict[str, Any],
    tool_name: str,
) -> dict[str, Any]:
    tool_object = api_response.get(tool_name)

    if isinstance(tool_object, dict):
        return tool_object

    tools_object = api_response.get("tools", {})

    if isinstance(tools_object, dict):
        nested_tool = tools_object.get(tool_name)

        if isinstance(nested_tool, dict):
            return nested_tool

    return {}


def get_rows(
    api_response: dict[str, Any],
    tool_name: str,
) -> list[dict[str, Any]]:
    results = api_response.get("results", {})

    if isinstance(results, dict):
        result_rows = results.get(tool_name)

        if isinstance(result_rows, list):
            return [
                row
                for row in result_rows
                if isinstance(row, dict)
            ]

    rows_object = api_response.get("rows", {})

    if isinstance(rows_object, dict):
        row_values = rows_object.get(tool_name)

        if isinstance(row_values, list):
            return [
                row
                for row in row_values
                if isinstance(row, dict)
            ]

    tool_object = _get_tool_object(
        api_response,
        tool_name,
    )

    fallback_rows = tool_object.get("rows", [])

    if isinstance(fallback_rows, list):
        return [
            row
            for row in fallback_rows
            if isinstance(row, dict)
        ]

    return []


def _normalize_tool_result(
    api_response: dict[str, Any],
    tool_name: str,
) -> dict[str, Any]:
    tool_object = _get_tool_object(
        api_response,
        tool_name,
    )

    rows = get_rows(
        api_response,
        tool_name,
    )

    error = tool_object.get("error")
    explicit_status = _normalize_status_value(
        tool_object.get("status")
    )

    source_status = _normalize_status_value(
        api_response.get("status")
    )

    if error:
        if explicit_status == "partial":
            status = "partial"
        else:
            status = "failed"

    elif explicit_status != "unknown":
        status = explicit_status

    elif source_status == "complete":
        status = (
            "complete"
            if rows
            else "complete_no_hits"
        )

    elif source_status == "failed":
        status = "failed"

    else:
        # A row count alone is not proof that the tool completed.
        # In an old or partial API response, preserve uncertainty
        # rather than converting an absent result into "no hits."
        status = "unknown"

    if status == "complete" and not rows:
        status = "complete_no_hits"

    row_count_reported = tool_object.get(
        "row_count_returned"
    )

    if not isinstance(row_count_reported, int):
        row_count_reported = None

    command = tool_object.get("command", [])

    if not isinstance(command, list):
        command = []

    return {
        "tool": tool_name,
        "status": status,
        "row_count": len(rows),
        "row_count_reported": row_count_reported,
        "rows": rows,
        "error": error,
        "command": command,
        "log_file": tool_object.get("log_file"),
        "output_file": tool_object.get("output_file"),
    }


def _aggregate_tool_status(
    tool_results: dict[str, dict[str, Any]],
) -> str:
    statuses = [
        result.get("status", "unknown")
        for result in tool_results.values()
    ]

    if not statuses:
        return "unknown"

    if all(
        status in SUCCESSFUL_TOOL_STATUSES
        for status in statuses
    ):
        return "complete"

    successful_count = sum(
        status in SUCCESSFUL_TOOL_STATUSES
        for status in statuses
    )

    if successful_count:
        return "partial"

    if any(status == "partial" for status in statuses):
        return "partial"

    if any(status == "failed" for status in statuses):
        return "failed"

    if all(status == "skipped" for status in statuses):
        return "skipped"

    return "unknown"


def parse_cpu_response(
    api_response: dict[str, Any],
) -> dict[str, Any]:
    tool_results = {
        tool_name: _normalize_tool_result(
            api_response,
            tool_name,
        )
        for tool_name in TOOL_NAMES
    }

    aggregate_status = _aggregate_tool_status(
        tool_results
    )

    result_counts = {
        tool_name: tool_results[tool_name]["row_count"]
        for tool_name in TOOL_NAMES
    }

    rows = {
        tool_name: tool_results[tool_name]["rows"]
        for tool_name in TOOL_NAMES
    }

    tool_statuses = {
        tool_name: tool_results[tool_name]["status"]
        for tool_name in TOOL_NAMES
    }

    tool_errors = {
        tool_name: tool_results[tool_name]["error"]
        for tool_name in TOOL_NAMES
        if tool_results[tool_name].get("error")
    }

    return {
        "job_id": api_response.get("job_id"),
        "status": aggregate_status,
        "source_status": api_response.get("status"),
        "files": api_response.get("files", {}),
        "result_counts": result_counts,
        "tool_statuses": tool_statuses,
        "tool_errors": tool_errors,
        "tools": tool_results,
        "rows": rows,
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

    raw_output = (
        output_dir
        / "homolog_search_raw.json"
    )
    summary_output = (
        output_dir
        / "homolog_search_summary.json"
    )
    error_output = (
        output_dir
        / "homolog_search_error.json"
    )

    try:
        raw_response = post_fasta(
            api_url=api_url,
            fasta_path=fasta_path,
            timeout_seconds=timeout_seconds,
        )

        parsed = parse_cpu_response(
            raw_response
        )

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
            source_fasta=str(
                fasta_path.resolve()
            ),
        )

        target_outputs = (
            write_target_evidence_outputs(
                target_evidence,
                output_dir,
            )
        )

        print(
            "[CPU] Annotation status: "
            f"{parsed.get('status')}",
            flush=True,
        )

        for tool_name in TOOL_NAMES:
            tool_result = parsed["tools"][tool_name]

            print(
                f"[CPU] {tool_name.upper()}: "
                f"status={tool_result['status']}; "
                f"rows={tool_result['row_count']}",
                flush=True,
            )

        return {
            # "ok" means the API request and response parsing
            # completed. It does not mean every annotation tool
            # completed successfully.
            "status": "ok",
            "cpu_status": parsed.get("status"),
            "api_url": api_url,
            "fasta_path": str(
                fasta_path.resolve()
            ),
            "raw_output": str(raw_output),
            "summary_output": str(summary_output),
            "target_evidence": target_outputs.get(
                "target_evidence"
            ),
            "target_evidence_report": (
                target_outputs.get(
                    "target_evidence_report"
                )
            ),
            "result_counts": parsed.get(
                "result_counts",
                {},
            ),
            "tool_statuses": parsed.get(
                "tool_statuses",
                {},
            ),
            "tool_errors": parsed.get(
                "tool_errors",
                {},
            ),
        }

    except Exception as error:
        payload = {
            "status": "error",
            "api_url": api_url,
            "fasta_path": str(
                fasta_path.resolve()
            ),
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
            "cpu_status": "failed",
            "api_url": api_url,
            "fasta_path": str(
                fasta_path.resolve()
            ),
            "error_output": str(error_output),
            "error_type": type(error).__name__,
            "error": str(error),
            "result_counts": {},
            "tool_statuses": {},
            "tool_errors": {},
        }
