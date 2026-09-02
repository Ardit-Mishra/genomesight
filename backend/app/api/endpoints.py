from fastapi import APIRouter, HTTPException, UploadFile, File, Form
from typing import List, Optional
from app.api.schemas import (
    AnalyzeRequest,
    AnalyzeResponse,
    AlignRequest,
    AlignResponse,
    TranslateRequest,
    TranslateResponse,
    CodonRequest,
    CodonResponse,
)
from app.core.io_fasta import parse_fasta_string
from app.core.alignment import perform_alignment
from app.core.codons import calculate_codon_usage, calculate_rscu
from app.core.kmer_native import count_kmers
from app.core.pinger import self_ping_status
from Bio.Seq import Seq

router = APIRouter(prefix="/api", tags=["GenomeSight Analysis"])

@router.post("/analyze", response_model=AnalyzeResponse)
async def analyze_sequence(payload: AnalyzeRequest):
    content = ""
    if payload.sequence:
        content = payload.sequence
    elif payload.file_content:
        content = payload.file_content
    else:
        raise HTTPException(status_code=400, detail="Either sequence or file_content must be provided.")

    records = parse_fasta_string(content)
    if not records:
        raise HTTPException(status_code=400, detail="No valid FASTA/FASTQ records found in input.")

    # Aggregate statistics across records
    total_length = sum(len(rec.seq) for rec in records)
    first_seq = str(records[0].seq).upper()

    # Calculate GC content
    g_count = first_seq.count('G')
    c_count = first_seq.count('C')
    gc_content = ((g_count + c_count) / len(first_seq) * 100.0) if len(first_seq) > 0 else 0.0

    # Nucleotide composition
    composition = {
        "A": first_seq.count('A'),
        "T": first_seq.count('T'),
        "G": g_count,
        "C": c_count,
        "N": first_seq.count('N')
    }

    statistics = {
        "record_count": len(records),
        "total_length": total_length,
        "primary_id": records[0].id,
        "primary_description": records[0].description,
        "gc_content_percentage": round(gc_content, 2),
        "composition": composition,
        "molecular_weight_approx": round(total_length * 330.0, 2) # rough estimate
    }

    # Quality statistics if FASTQ format
    quality_stats = None
    if any(hasattr(rec, "letter_annotations") and "phred_quality" in rec.letter_annotations for rec in records):
        qualities = [q for rec in records for q in rec.letter_annotations.get("phred_quality", [])]
        if qualities:
            quality_stats = {
                "mean_quality_score": round(sum(qualities) / len(qualities), 2),
                "max_quality": max(qualities),
                "min_quality": min(qualities)
            }

    # K-mer analysis
    seq_list = [str(rec.seq) for rec in records]
    kmer_counts, engine = count_kmers(seq_list, k=3)
    kmer_analysis = {
        "k": 3,
        "engine": engine,
        "counts": dict(list(kmer_counts.items())[:50]) # top 50 for performance
    }

    return AnalyzeResponse(
        success=True,
        statistics=statistics,
        quality_statistics=quality_stats,
        kmer_analysis=kmer_analysis
    )

@router.post("/align", response_model=AlignResponse)
async def align_sequences(payload: AlignRequest):
    if not payload.seq1 or not payload.seq2:
        raise HTTPException(status_code=400, detail="Both seq1 and seq2 must be provided.")

    result = perform_alignment(payload.seq1, payload.seq2, mode=payload.mode)
    return AlignResponse(
        success=True,
        score=result["score"],
        aligned_seq1=result["aligned_seq1"],
        aligned_seq2=result["aligned_seq2"],
        match_line=result["match_line"],
        identity_percentage=result["identity_percentage"]
    )

@router.post("/translate", response_model=TranslateResponse)
async def translate_sequence(payload: TranslateRequest):
    if not payload.sequence:
        raise HTTPException(status_code=400, detail="Sequence must be provided.")

    if payload.frame not in [1, 2, 3]:
        raise HTTPException(status_code=400, detail="Reading frame must be 1, 2, or 3.")

    seq_str = payload.sequence.upper().replace("U", "T")
    # Adjust for reading frame
    offset = payload.frame - 1
    trimmed_seq = seq_str[offset:]

    try:
        biopython_seq = Seq(trimmed_seq)
        protein = str(biopython_seq.translate(to_stop=False))
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Translation failed: {str(e)}")

    return TranslateResponse(
        success=True,
        protein_sequence=protein,
        frame=payload.frame
    )

@router.post("/codons", response_model=CodonResponse)
async def codon_usage(payload: CodonRequest):
    if not payload.sequence:
        raise HTTPException(status_code=400, detail="Sequence must be provided.")

    counts = calculate_codon_usage(payload.sequence)
    rscu = calculate_rscu(payload.sequence)

    return CodonResponse(
        success=True,
        codon_counts=counts,
        rscu=rscu
    )

@router.get("/health")
async def health_check():
    return {
        "status": self_ping_status(),
        "service": "genomesight-backend",
        "version": "2.0.0"
    }
