from Bio import Align
from typing import Any, Dict


def perform_alignment(seq1: str, seq2: str, mode: str = "global") -> Dict[str, Any]:
    """Perform pairwise alignment and return clean, machine-readable rows."""
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return {
            "score": 0.0,
            "aligned_seq1": seq1,
            "aligned_seq2": seq2,
            "match_line": "",
            "identity_percentage": 0.0,
        }

    best_alignment = alignments[0]
    aligned_seq1, aligned_seq2 = _gapped_sequences(best_alignment, seq1, seq2)
    match_line = "".join(
        "|" if base1 == base2 and base1 != "-" else " "
        for base1, base2 in zip(aligned_seq1, aligned_seq2)
    )
    identity_percentage = (match_line.count("|") / max(len(aligned_seq1), 1)) * 100.0

    return {
        "score": float(best_alignment.score),
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
        "match_line": match_line,
        "identity_percentage": float(identity_percentage),
    }


def _gapped_sequences(alignment: Align.Alignment, seq1: str, seq2: str) -> tuple[str, str]:
    """Build gapped rows from PairwiseAligner's coordinate path."""
    pieces1: list[str] = []
    pieces2: list[str] = []
    coordinates = alignment.coordinates

    for index in range(coordinates.shape[1] - 1):
        start1, end1 = coordinates[0, index], coordinates[0, index + 1]
        start2, end2 = coordinates[1, index], coordinates[1, index + 1]
        segment1 = seq1[start1:end1]
        segment2 = seq2[start2:end2]

        if len(segment1) < len(segment2):
            segment1 = "-" * (len(segment2) - len(segment1)) + segment1
        elif len(segment2) < len(segment1):
            segment2 = "-" * (len(segment1) - len(segment2)) + segment2

        pieces1.append(segment1)
        pieces2.append(segment2)

    return "".join(pieces1), "".join(pieces2)
