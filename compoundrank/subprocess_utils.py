from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Mapping, Sequence


def resolve_executable(value: str, label: str | None = None) -> str:
    expanded = str(Path(value).expanduser())
    resolved = shutil.which(expanded)
    if resolved:
        return resolved
    candidate = Path(expanded)
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate.resolve())
    raise RuntimeError(f"{label or value} executable was not found: {value}")


def run_command(
    command: Sequence[str],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    capture_output: bool = True,
) -> subprocess.CompletedProcess[str]:
    print("$", " ".join(shlex.quote(str(item)) for item in command), flush=True)
    completed = subprocess.run(
        [str(item) for item in command],
        cwd=str(cwd) if cwd else None,
        env=dict(env) if env else None,
        text=True,
        capture_output=capture_output,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code "
            f"{completed.returncode}: {' '.join(command)}\n\n"
            f"STDOUT:\n{completed.stdout or ''}\n\n"
            f"STDERR:\n{completed.stderr or ''}"
        )
    return completed
