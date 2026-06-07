from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import JSONResponse

from cpu_server.runners import run_cpu_analysis
from cpu_server.schemas import FastaAnalysisResponse

from .settings import (
    DATA_ROOT,
    DEPLOY_ROOT,
    INTERPRO_DATA_ROOT,
    JOBS_ROOT,
    REPO_ROOT,
    VOGDB_DATA_ROOT,
)


app = FastAPI(
    title="CompoundRank CPU Server",
    description="CPU-only API for FASTA parsing, CDD, InterPro, and VOGDB analysis.",
    version="0.1.0",
)

@app.get("/debug/paths")
def debug_paths():
    return {
        "repo_root": str(REPO_ROOT),
        "deploy_root": str(DEPLOY_ROOT),
        "data_root": str(DATA_ROOT),
        "jobs_root": str(JOBS_ROOT),
        "interpro_data_root": str(INTERPRO_DATA_ROOT),
        "vogdb_data_root": str(VOGDB_DATA_ROOT),
        "interpro_data_exists": INTERPRO_DATA_ROOT.exists(),
        "vogdb_data_exists": VOGDB_DATA_ROOT.exists(),
        "jobs_root_exists": JOBS_ROOT.exists(),
    }

@app.get("/health")
def health():
    return {
        "status": "ok",
        "service": "compoundrank-cpu-server",
        "cpu_steps": ["fasta_parse", "cdd", "interpro", "vogdb"],
        "gpu_steps_excluded": ["colabfold", "gnina"],
    }


@app.post("/analyze/fasta", response_model=FastaAnalysisResponse)
def analyze_fasta(file: UploadFile = File(...)):
    filename = file.filename or ""

    if not filename.endswith((".fa", ".faa", ".fasta", ".txt")):
        raise HTTPException(
            status_code=400,
            detail="Upload must be a FASTA-like file: .fa, .faa, .fasta, or .txt",
        )

    try:
        result = run_cpu_analysis(file)
        return JSONResponse(content=result)

    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"CPU analysis failed: {exc}")