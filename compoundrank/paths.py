from __future__ import annotations

import hashlib
from pathlib import Path


def repository_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def require_absolute_external_file(value: str | Path, label: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise ValueError(f"{label} must be an absolute path: {path}")
    path = path.resolve()
    if _is_relative_to(path, repository_root()):
        raise ValueError(
            f"{label} must be outside the repository. Received: {path}"
        )
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if path.stat().st_size == 0:
        raise RuntimeError(f"{label} is empty: {path}")
    return path


def require_absolute_external_dir(
    value: str | Path,
    label: str,
    create: bool = False,
) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise ValueError(f"{label} must be an absolute path: {path}")
    path = path.resolve()
    if _is_relative_to(path, repository_root()):
        raise ValueError(
            f"{label} must be outside the repository. Received: {path}"
        )
    if create:
        path.mkdir(parents=True, exist_ok=True)
    if not path.is_dir():
        raise NotADirectoryError(f"{label} is not a directory: {path}")
    return path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def content_cache_key(*parts: str | Path) -> str:
    digest = hashlib.sha256()
    for part in parts:
        if isinstance(part, Path):
            digest.update(file_sha256(part).encode("ascii"))
        else:
            digest.update(str(part).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()[:24]


def sanitize_name(value: str) -> str:
    cleaned = []
    previous_underscore = False
    for character in value.strip():
        if character.isalnum() or character in {"-", "."}:
            cleaned.append(character.lower())
            previous_underscore = False
        else:
            if not previous_underscore:
                cleaned.append("_")
                previous_underscore = True
    result = "".join(cleaned).strip("_.-")
    return result or "ligand"
