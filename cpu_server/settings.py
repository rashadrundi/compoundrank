import os
from pathlib import Path


# Expected droplet layout:
#
# /opt/compoundrank/
#   data/
#     interpro_data/
#     vogdb_data/
#   jobs/
#   repo/
#     compoundrank/
#       cpu_server/
#       av_pipeline/
#       ...


REPO_ROOT = Path(
    os.getenv(
        "COMPOUNDRANK_REPO_ROOT",
        Path(__file__).resolve().parents[1],
    )
).resolve()


# If repo is /opt/compoundrank/repo/compoundrank,
# then REPO_ROOT.parent is /opt/compoundrank/repo
# and REPO_ROOT.parents[1] is /opt/compoundrank.
if REPO_ROOT.parent.name == "repo":
    DEFAULT_DEPLOY_ROOT = REPO_ROOT.parents[1]
else:
    DEFAULT_DEPLOY_ROOT = REPO_ROOT


DEPLOY_ROOT = Path(
    os.getenv(
        "COMPOUNDRANK_DEPLOY_ROOT",
        DEFAULT_DEPLOY_ROOT,
    )
).resolve()


DATA_ROOT = Path(
    os.getenv(
        "COMPOUNDRANK_DATA_ROOT",
        DEPLOY_ROOT / "data",
    )
).resolve()


JOBS_ROOT = Path(
    os.getenv(
        "COMPOUNDRANK_JOBS_ROOT",
        DEPLOY_ROOT / "jobs",
    )
).resolve()


INTERPRO_DATA_ROOT = Path(
    os.getenv(
        "INTERPRO_DATA_ROOT",
        DATA_ROOT / "interpro_data",
    )
).resolve()


VOGDB_DATA_ROOT = Path(
    os.getenv(
        "VOGDB_DATA_ROOT",
        DATA_ROOT / "vogdb_data",
    )
).resolve()


def ensure_runtime_dirs() -> None:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    JOBS_ROOT.mkdir(parents=True, exist_ok=True)


def runtime_env() -> dict[str, str]:
    return {
        "COMPOUNDRANK_REPO_ROOT": str(REPO_ROOT),
        "COMPOUNDRANK_DEPLOY_ROOT": str(DEPLOY_ROOT),
        "COMPOUNDRANK_DATA_ROOT": str(DATA_ROOT),
        "COMPOUNDRANK_JOBS_ROOT": str(JOBS_ROOT),
        "INTERPRO_DATA_ROOT": str(INTERPRO_DATA_ROOT),
        "VOGDB_DATA_ROOT": str(VOGDB_DATA_ROOT),
    }


def apply_runtime_env() -> None:
    os.environ.update(runtime_env())