from Bio import Align
from typing import Dict, Any

def perform_alignment(seq1: str, seq2: str, mode: str = "global") -> Dict[str, Any]:
    """Perform pairwise sequence alignment using Bio.Align."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global" if mode.lower() == "global" else "local"
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return {
            "score": 0.0,
            "aligned_seq1": seq1,
            "aligned_seq2": seq2,
            "match_line": "",
            "identity_percentage": 0.0
        }

    best_alignment = alignments[0]
    lines = format(best_alignment).split("\n")

    aligned_s1 = lines[0] if len(lines) > 0 else seq1
    match_ln = lines[1] if len(lines) > 1 else ""
    aligned_s2 = lines[2] if len(lines) > 2 else seq2

    matches = match_ln.count("|")
    total_len = len(match_ln) if len(match_ln) > 0 else 1
    identity = (matches / total_len) * 100.0

    return {
        "score": float(best_alignment.score),
        "aligned_seq1": aligned_s1,
        "aligned_seq2": aligned_s2,
        "match_line": match_ln,
        "identity_percentage": float(identity)
    }
