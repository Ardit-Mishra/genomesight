from pydantic import BaseModel, Field
from typing import Dict, List, Any, Optional

class AnalyzeRequest(BaseModel):
    sequence: Optional[str] = Field(None, description="Raw DNA/RNA sequence string")
    file_content: Optional[str] = Field(None, description="FASTA/FASTQ file content as string")

class AnalyzeResponse(BaseModel):
    success: bool
    statistics: Dict[str, Any]
    quality_statistics: Optional[Dict[str, Any]] = None
    kmer_analysis: Dict[str, Any]

class AlignRequest(BaseModel):
    seq1: str = Field(..., description="First sequence")
    seq2: str = Field(..., description="Second sequence")
    mode: str = Field("global", description="Alignment mode: 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman)")

class AlignResponse(BaseModel):
    success: bool
    score: float
    aligned_seq1: str
    aligned_seq2: str
    match_line: str
    identity_percentage: float

class TranslateRequest(BaseModel):
    sequence: str = Field(..., description="Nucleotide sequence")
    frame: int = Field(1, description="Reading frame (1, 2, or 3)")

class TranslateResponse(BaseModel):
    success: bool
    protein_sequence: str
    frame: int

class CodonRequest(BaseModel):
    sequence: str = Field(..., description="Coding sequence")

class CodonResponse(BaseModel):
    success: bool
    codon_counts: Dict[str, int]
    rscu: Dict[str, float]
