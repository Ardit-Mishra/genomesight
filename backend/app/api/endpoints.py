from fastapi import APIRouter, HTTPException
from Bio.Seq import Seq

from app.api.schemas import (
    AlignRequest,
    AlignResponse,
    AnalyzeRequest,
    AnalyzeResponse,
    CodonRequest,
    CodonResponse,
    EnzymeListResponse,
    MotifRequest,
    MotifResponse,
    OrfRequest,
    OrfResponse,
    RestrictionRequest,
    TranslateRequest,
    TranslateResponse,
)
from app.core.alignment import perform_alignment
from app.core.build_info import build_info
from app.core.codons import calculate_codon_usage, calculate_rscu
from app.core.io_fasta import parse_fasta_string
from app.core.kmer_native import backend_status as kmer_backend_status, count_kmers
from app.core.motifs import (
    get_enzyme_list,
    iupac_to_regex,
    search_motif,
    search_restriction_site,
)
from app.core.orf import find_orfs, get_orf_summary
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
        raise HTTPException(status_code=400, detail="No valid FASTA, FASTQ or GenBank records found in input.")

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
    """Liveness plus the deployed commit.

    The commit is the point: without it, telling a stale deployment from a code
    defect means comparing behaviour endpoint by endpoint. With it, it is one
    field against `git rev-parse HEAD`.
    """
    info = build_info()
    return {
        "status": self_ping_status(),
        "service": "genomesight-backend",
        "version": "2.0.0",
        "commit": info["commit"],
        "commit_short": info["commit_short"],
        "branch": info["branch"],
        "commit_source": info["source"],
        "kmer_engine": kmer_backend_status()["backend"],
    }


def _records_from(payload_sequence: str):
    """Parse input to SeqRecords, accepting either FASTA/FASTQ or a bare sequence.

    Shared by the ORF and motif endpoints so all three accept the same inputs the
    analyze endpoint does, rather than each inventing its own rules.
    """
    content = payload_sequence or ""
    if not content.strip():
        raise HTTPException(status_code=400, detail="A sequence must be provided.")
    if not content.lstrip().startswith((">", "@")):
        try:
            content = wrap_raw_as_fasta(content)
        except ValueError as error:
            raise _http_422(error) from error
    records = parse_fasta_string(content)
    if not records:
        raise HTTPException(status_code=400, detail="No valid FASTA, FASTQ or GenBank records found in input.")
    return records


@router.post("/orfs", response_model=OrfResponse)
def detect_orfs(payload: OrfRequest):
    """Six-frame open reading frame detection."""
    records = _records_from(payload.sequence)
    try:
        sequence = "".join(normalize_sequence(str(record.seq)) for record in records)
    except ValueError as error:
        raise _http_422(error) from error

    orfs = find_orfs(
        sequence,
        min_length=payload.min_length,
        include_reverse=payload.include_reverse,
        use_alternative_starts=payload.use_alternative_starts,
    )
    summary = get_orf_summary(orfs)
    # The summary carries the ORF dataclass itself under 'longest_orf', which is not
    # JSON-serialisable; the full record is already in the orfs list.
    summary.pop("longest_orf", None)

    return {
        "success": True,
        "coordinate_note": (
            "Minus-strand coordinates are given against the reverse-complement "
            "sequence, not the original strand."
        ),
        "total": len(orfs),
        "summary": summary,
        "orfs": [
            {
                "start": orf.start,
                "end": orf.end,
                "length_nt": orf.length_nt,
                "length_aa": orf.length_aa,
                "frame": orf.frame,
                "strand": orf.strand,
                "start_codon": orf.start_codon,
                "stop_codon": orf.stop_codon,
                "protein": orf.protein,
                "gc_content": round(orf.gc_content, 2),
            }
            for orf in orfs
        ],
    }


@router.post("/motifs", response_model=MotifResponse)
def find_motifs(payload: MotifRequest):
    """Search for an IUPAC nucleotide pattern (R, Y, N, ... are expanded)."""
    records = _records_from(payload.sequence)
    matches = search_motif(records, payload.pattern, context_length=payload.context_length)
    return {
        "success": True,
        "pattern": payload.pattern.upper(),
        # Returning the compiled regex makes the ambiguity expansion inspectable
        # instead of asking the user to trust it.
        "regex": iupac_to_regex(payload.pattern),
        "total": len(matches),
        "matches": [
            {
                "pattern": m.pattern,
                "sequence_id": m.sequence_id,
                "start_1based": m.start_1based,
                "end": m.end,
                "matched_sequence": m.matched_sequence,
                "context": m.context,
            }
            for m in matches
        ],
    }


@router.post("/restriction-sites", response_model=MotifResponse)
def find_restriction_sites(payload: RestrictionRequest):
    """Search for a named restriction enzyme's recognition site."""
    records = _records_from(payload.sequence)
    try:
        matches = search_restriction_site(records, payload.enzyme)
    except ValueError as error:
        raise _http_422(error) from error

    site = get_enzyme_list()[payload.enzyme]
    return {
        "success": True,
        "pattern": f"{payload.enzyme} ({site})",
        "regex": iupac_to_regex(site),
        "total": len(matches),
        "matches": [
            {
                "pattern": m.pattern,
                "sequence_id": m.sequence_id,
                "start_1based": m.start_1based,
                "end": m.end,
                "matched_sequence": m.matched_sequence,
                "context": m.context,
            }
            for m in matches
        ],
    }


@router.get("/enzymes", response_model=EnzymeListResponse)
def list_enzymes():
    """Recognition sites for the restriction enzymes this service knows."""
    return {"success": True, "enzymes": get_enzyme_list()}
