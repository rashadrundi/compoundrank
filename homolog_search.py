from __future__ import annotations

import json
import uuid
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


DEFAULT_API_URL = "http://161.35.0.191:8000/analyze/fasta"


def build_multipart_body(file_path: Path, field_name: str = "file") -> tuple[bytes, str]:
    boundary = f"----CompoundRankBoundary{uuid.uuid4().hex}"
    filename = file_path.name
    file_bytes = file_path.read_bytes()

    body = b"".join([
        f"--{boundary}\r\n".encode(),
        (
            f'Content-Disposition: form-data; name="{field_name}"; '
            f'filename="{filename}"\r\n'
        ).encode(),
        b"Content-Type: application/octet-stream\r\n\r\n",
        file_bytes,
        b"\r\n",
        f"--{boundary}--\r\n".encode(),
    ])

    content_type = f"multipart/form-data; boundary={boundary}"
    return body, content_type


def post_fasta(api_url: str, fasta_path: Path, timeout_seconds: int = 3600) -> dict[str, Any]:
    body, content_type = build_multipart_body(fasta_path)

    request = urllib.request.Request(
        api_url,
        data=body,
        method="POST",
        headers={
            "Content-Type": content_type,
            "Content-Length": str(len(body)),
        },
    )

    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            response_text = response.read().decode("utf-8", errors="replace")
            status_code = response.status

    except urllib.error.HTTPError as error:
        error_body = error.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP error from CPU API: {error.code}\n{error_body}") from error

    except urllib.error.URLError as error:
        raise RuntimeError(f"Could not reach CPU API: {error}") from error

    try:
        return json.loads(response_text)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"CPU API returned HTTP {status_code}, but response was not JSON:\n"
            f"{response_text[:1000]}"
        ) from error


def get_rows(api_response: dict[str, Any], tool_name: str) -> list[dict[str, Any]]:
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
