from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Mapping, Sequence


class CommandTimeoutError(RuntimeError):
    """Raised when an external command exceeds its time limit."""

    def __init__(
        self,
        *,
        command: Sequence[str],
        timeout_seconds: float,
        stdout: str = "",
        stderr: str = "",
    ) -> None:
        self.command = [str(item) for item in command]
        self.timeout_seconds = float(timeout_seconds)
        self.stdout = stdout
        self.stderr = stderr

        super().__init__(
            "Command timed out after "
            f"{self.timeout_seconds:g} seconds: "
            f"{' '.join(self.command)}\n\n"
            f"STDOUT:\n{stdout}\n\n"
            f"STDERR:\n{stderr}"
        )


def _timeout_output(value: object) -> str:
    if value is None:
        return ""

    if isinstance(value, bytes):
        return value.decode(
            "utf-8",
            errors="replace",
        )

    return str(value)


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
    timeout_seconds: float | None = None,
) -> subprocess.CompletedProcess[str]:
    print(
        "$",
        " ".join(
            shlex.quote(str(item))
            for item in command
        ),
        flush=True,
    )

    try:
        completed = subprocess.run(
            [str(item) for item in command],
            cwd=str(cwd) if cwd else None,
            env=dict(env) if env else None,
            text=True,
            capture_output=capture_output,
            check=False,
            timeout=timeout_seconds,
        )

    except subprocess.TimeoutExpired as error:
        raise CommandTimeoutError(
            command=command,
            timeout_seconds=(
                timeout_seconds
                if timeout_seconds is not None
                else 0
            ),
            stdout=_timeout_output(error.stdout),
            stderr=_timeout_output(error.stderr),
        ) from error
    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code "
            f"{completed.returncode}: {' '.join(command)}\n\n"
            f"STDOUT:\n{completed.stdout or ''}\n\n"
            f"STDERR:\n{completed.stderr or ''}"
        )
    return completed
