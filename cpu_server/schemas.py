from pydantic import BaseModel
from typing import Any, Dict, List, Optional


class ToolResult(BaseModel):
    status: str
    output_file: Optional[str] = None
    rows: List[Dict[str, Any]] = []
    error: Optional[str] = None


class FastaAnalysisResponse(BaseModel):
    job_id: str
    status: str
    fasta_file: str
    cdd: ToolResult
    interpro: ToolResult
    vogdb: ToolResult