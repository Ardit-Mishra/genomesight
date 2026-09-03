from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field


class AnalyzeRequest(BaseModel):
    sequence: Optional[str] = Field(None, description="Raw DNA/RNA sequence string or FASTA/FASTQ content")
    file_content: Optional[str] = Field(None, description="FASTA/FASTQ file content as string")


class AnalyzeResponse(BaseModel):
    success: bool
    statistics: Dict[str, Any]
    quality_statistics: Optional[Dict[str, Any]] = None
    kmer_analysis: Dict[str, Any]


class AlignRequest(BaseModel):
    seq1: str = Field(..., description="First sequence")
    seq2: str = Field(..., description="Second sequence")
    mode: Literal["global", "local"] = Field("global", description="Needleman-Wunsch global or Smith-Waterman local alignment")


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


# ---------------------------------------------------------------------------
# ORF detection
# ---------------------------------------------------------------------------
class OrfRequest(BaseModel):
    sequence: str
    min_length: int = Field(default=100, ge=3, le=100_000)
    include_reverse: bool = True
    use_alternative_starts: bool = False


class OrfRecord(BaseModel):
    start: int
    end: int
    length_nt: int
    length_aa: int
    frame: int
    strand: str
    start_codon: str
    stop_codon: str
    protein: str
    gc_content: float


class OrfResponse(BaseModel):
    success: bool
    # Coordinates on the minus strand are reported against the reverse-complement
    # sequence, not the original. Stated here because a coordinate whose frame of
    # reference is unstated is the kind of number that gets silently misread.
    coordinate_note: str
    total: int
    summary: dict
    orfs: List[OrfRecord]


# ---------------------------------------------------------------------------
# Motif and restriction-site search
# ---------------------------------------------------------------------------
class MotifRequest(BaseModel):
    sequence: str
    pattern: str = Field(min_length=1, max_length=64)
    context_length: int = Field(default=10, ge=0, le=100)


class RestrictionRequest(BaseModel):
    sequence: str
    enzyme: str


class MotifMatchRecord(BaseModel):
    pattern: str
    sequence_id: str
    start_1based: int
    end: int
    matched_sequence: str
    context: str


class MotifResponse(BaseModel):
    success: bool
    pattern: str
    regex: str
    total: int
    matches: List[MotifMatchRecord]


class EnzymeListResponse(BaseModel):
    success: bool
    enzymes: dict
