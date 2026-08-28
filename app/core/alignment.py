"""
Pairwise Sequence Alignment Module

This module provides functions for aligning two sequences using global
(Needleman-Wunsch) and local (Smith-Waterman) alignment algorithms.
It wraps Biopython's PairwiseAligner for robust implementation.

Design Notes:
    - Uses Biopython PairwiseAligner for algorithm implementation
    - Supports both global and local alignment modes
    - Calculates identity percentage and generates alignment visualization
    - Returns detailed alignment statistics

Alignment Algorithms:
    Global: Aligns entire sequences end-to-end (Needleman-Wunsch)
    Local: Finds best matching subsequence (Smith-Waterman)
"""

from Bio import Align
from Bio.Seq import Seq
from typing import Dict, Any, Optional, Tuple, List
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class AlignmentResult:
    """Data class for alignment results."""
    score: float
    aligned_seq1: str
    aligned_seq2: str
    match_line: str
    identity: float
    identity_count: int
    alignment_length: int
    gaps: int
    mode: str


def align_sequences(
    seq1: str,
    seq2: str,
    mode: str = "global",
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    gap_open: float = -10.0,
    gap_extend: float = -0.5
) -> AlignmentResult:
    """
    Perform pairwise sequence alignment.
    
    Args:
        seq1: First sequence string
        seq2: Second sequence string
        mode: "global" or "local" alignment
        match_score: Score for matching bases
        mismatch_score: Penalty for mismatches
        gap_open: Penalty for opening a gap
        gap_extend: Penalty for extending a gap
    
    Returns:
        AlignmentResult with alignment details
    
    Example:
        >>> result = align_sequences("ATGC", "ATCC")
        >>> print(f"Identity: {result.identity:.1f}%")
    """
    seq1 = seq1.upper().replace(' ', '')
    seq2 = seq2.upper().replace(' ', '')
    
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    
    try:
        alignments = aligner.align(seq1, seq2)
        
        if not alignments:
            return _empty_result(mode)
        
        best_alignment = alignments[0]
        score = best_alignment.score
        
        aligned_seq1, aligned_seq2, match_line = _format_alignment(best_alignment)
        
        identity_count = sum(
            1 for a, b in zip(aligned_seq1, aligned_seq2) 
            if a == b and a != '-'
        )
        alignment_length = len(aligned_seq1)
        gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
        
        identity = (identity_count / alignment_length * 100) if alignment_length > 0 else 0
        
        logger.info(f"Alignment complete: {identity:.1f}% identity")
        
        return AlignmentResult(
            score=float(score),
            aligned_seq1=aligned_seq1,
            aligned_seq2=aligned_seq2,
            match_line=match_line,
            identity=identity,
            identity_count=identity_count,
            alignment_length=alignment_length,
            gaps=gaps,
            mode=mode
        )
        
    except Exception as e:
        logger.error(f"Alignment failed: {e}")
        return _empty_result(mode)


def _format_alignment(alignment) -> Tuple[str, str, str]:
    """
    Extract the two gapped, aligned rows from a Biopython Alignment.

    Indexing the alignment (`alignment[0]`, `alignment[1]`) yields exactly the
    aligned sequences, gaps included, and nothing else.

    This previously parsed `alignment.format()` and took lines 0 and 2. From
    Biopython 1.80 that output carries coordinate columns —

        target            0 ATGCGTACGTTAGC 14
                          0 ||||||.||||||| 14
        query             0 ATGCGTTCGTTAGC 14

    — so the labels and offsets were consumed as sequence characters. Every
    downstream number was computed over that text: two 14 bp sequences differing
    at one position reported 81.1% identity over a 37-position alignment instead
    of 92.9% over 14.
    """
    try:
        aligned_seq1 = str(alignment[0])
        aligned_seq2 = str(alignment[1])
    except Exception:  # pragma: no cover - defensive, older/odd Biopython
        logger.warning("Could not index alignment rows; falling back to ungapped sequences")
        try:
            aligned_seq1 = str(alignment.target)
            aligned_seq2 = str(alignment.query)
        except Exception:
            aligned_seq1 = ""
            aligned_seq2 = ""


    match_line = []
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == b and a != '-':
            match_line.append('|')
        elif a == '-' or b == '-':
            match_line.append(' ')
        else:
            match_line.append('.')
    
    return aligned_seq1, aligned_seq2, ''.join(match_line)


def _empty_result(mode: str) -> AlignmentResult:
    """Return empty result when alignment fails."""
    return AlignmentResult(
        score=0.0,
        aligned_seq1="",
        aligned_seq2="",
        match_line="",
        identity=0.0,
        identity_count=0,
        alignment_length=0,
        gaps=0,
        mode=mode
    )


def calculate_identity(seq1: str, seq2: str) -> float:
    """
    Calculate simple identity percentage between two sequences.
    
    Compares sequences position by position without gaps.
    Sequences should be pre-aligned or same length.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
    
    Returns:
        Identity percentage (0-100)
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / min_len) * 100


def format_alignment_display(
    result: AlignmentResult,
    line_width: int = 60
) -> str:
    """
    Format alignment result for text display.
    
    Args:
        result: AlignmentResult object
        line_width: Characters per line
    
    Returns:
        Formatted multi-line alignment string
    """
    lines = []
    lines.append(f"Alignment Mode: {result.mode.capitalize()}")
    lines.append(f"Score: {result.score:.1f}")
    lines.append(f"Identity: {result.identity:.1f}% ({result.identity_count}/{result.alignment_length})")
    lines.append(f"Gaps: {result.gaps}")
    lines.append("")
    
    seq1 = result.aligned_seq1
    seq2 = result.aligned_seq2
    match = result.match_line
    
    for i in range(0, len(seq1), line_width):
        lines.append(f"Query:  {seq1[i:i+line_width]}")
        lines.append(f"        {match[i:i+line_width]}")
        lines.append(f"Subject:{seq2[i:i+line_width]}")
        lines.append("")
    
    return '\n'.join(lines)
