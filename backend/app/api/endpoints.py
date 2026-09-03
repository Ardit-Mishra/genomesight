from fastapi import APIRouter, HTTPException
from Bio.Seq import Seq

from app.api.schemas import (
    AlignRequest,
    AlignResponse,
    AnalyzeRequest,
    AnalyzeResponse,
    CodonRequest,
    CodonResponse,
    TranslateRequest,
    TranslateResponse,
)
from app.core.alignment import perform_alignment
from app.core.codons import calculate_codon_usage, calculate_rscu
from app.core.io_fasta import parse_fasta_string
from app.core.kmer_native import count_kmers
from app.core.pinger import self_ping_status
from app.core.validation import normalize_sequence, wrap_raw_as_fasta

router = APIRouter(prefix="/api", tags=["GenomeSight Analysis"])


def _http_422(error: ValueError) -> HTTPException:
    return HTTPException(status_code=422, detail=str(error))


@router.post("/analyze", response_model=AnalyzeResponse)
def analyze_sequence(payload: AnalyzeRequest):
    if not payload.sequence and not payload.file_content:
        raise HTTPException(status_code=400, detail="Either sequence or file_content must be provided.")

    content = payload.file_content or payload.sequence or ""
    if payload.sequence and not payload.file_content and not content.lstrip().startswith((">", "@")):
        try:
            content = wrap_raw_as_fasta(content)
        except ValueError as error:
            raise _http_422(error) from error

    records = parse_fasta_string(content)
    if not records:
        raise HTTPException(status_code=400, detail="No valid FASTA/FASTQ records found in input.")

    try:
        sequences = [normalize_sequence(str(record.seq)) for record in records]
    except ValueError as error:
        raise _http_422(error) from error

    total_length = sum(len(sequence) for sequence in sequences)
    combined_sequence = "".join(sequences)
    g_count = combined_sequence.count("G")
    c_count = combined_sequence.count("C")
    composition = {
        "A": combined_sequence.count("A"),
        "T": combined_sequence.count("T"),
        "G": g_count,
        "C": c_count,
        "N": combined_sequence.count("N"),
    }

    quality_stats = None
    qualities = [
        quality
        for record in records
        for quality in record.letter_annotations.get("phred_quality", [])
    ]
    if qualities:
        quality_stats = {
            "mean_quality_score": round(sum(qualities) / len(qualities), 2),
            "max_quality": max(qualities),
            "min_quality": min(qualities),
        }

    kmer_counts, engine = count_kmers(sequences, k=3)
    return AnalyzeResponse(
        success=True,
        statistics={
            "record_count": len(records),
            "total_length": total_length,
            "primary_id": records[0].id,
            "primary_description": records[0].description,
            "gc_content_percentage": round(((g_count + c_count) / total_length * 100) if total_length else 0, 2),
            "composition": composition,
            "molecular_weight_approx": round(total_length * 330.0, 2),
        },
        quality_statistics=quality_stats,
        kmer_analysis={
            "k": 3,
            "engine": engine,
            "counts": dict(sorted(kmer_counts.items(), key=lambda item: (-item[1], item[0]))[:50]),
        },
    )


@router.post("/align", response_model=AlignResponse)
def align_sequences(payload: AlignRequest):
    try:
        seq1 = normalize_sequence(payload.seq1)
        seq2 = normalize_sequence(payload.seq2)
    except ValueError as error:
        raise _http_422(error) from error
    return AlignResponse(success=True, **perform_alignment(seq1, seq2, mode=payload.mode))


@router.post("/translate", response_model=TranslateResponse)
def translate_sequence(payload: TranslateRequest):
    if payload.frame not in [1, 2, 3]:
        raise HTTPException(status_code=400, detail="Reading frame must be 1, 2, or 3.")
    try:
        sequence = normalize_sequence(payload.sequence)
    except ValueError as error:
        raise _http_422(error) from error
    protein = str(Seq(sequence[payload.frame - 1:]).translate(to_stop=False))
    return TranslateResponse(success=True, protein_sequence=protein, frame=payload.frame)


@router.post("/codons", response_model=CodonResponse)
def codon_usage(payload: CodonRequest):
    try:
        sequence = normalize_sequence(payload.sequence)
    except ValueError as error:
        raise _http_422(error) from error
    return CodonResponse(
        success=True,
        codon_counts=calculate_codon_usage(sequence),
        rscu=calculate_rscu(sequence),
    )


@router.get("/health")
def health_check():
    return {"status": self_ping_status(), "service": "genomesight-backend", "version": "2.0.0"}
