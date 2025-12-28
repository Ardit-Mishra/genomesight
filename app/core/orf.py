"""
Open Reading Frame (ORF) Detection Module

This module provides functions for finding and analyzing Open Reading Frames
in DNA sequences. ORFs are regions that potentially encode proteins.

Design Notes:
    - Scans all three forward reading frames
    - Optionally scans reverse complement frames
    - Supports custom start/stop codon sets
    - Uses standard genetic code for translation

ORF Detection Algorithm:
    1. Scan sequence in triplets (codons) starting at position 0, 1, or 2
    2. Find start codons (default: ATG)
    3. Continue until stop codon (TAA, TAG, TGA)
    4. Report ORFs meeting minimum length threshold
"""

from typing import List, Dict, Any, Optional, Set
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


STANDARD_START_CODONS = {'ATG'}
ALTERNATIVE_START_CODONS = {'GTG', 'TTG', 'CTG'}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


@dataclass
class ORF:
    """Data class representing an Open Reading Frame."""
    start: int
    end: int
    length_nt: int
    length_aa: int
    frame: int
    strand: str
    start_codon: str
    stop_codon: str
    sequence: str
    protein: str
    gc_content: float


def translate_sequence(sequence: str, frame: int = 0) -> str:
    """
    Translate a DNA sequence to protein.
    
    Uses the standard genetic code. Translation stops at the first
    stop codon encountered.
    
    Args:
        sequence: DNA sequence string
        frame: Reading frame offset (0, 1, or 2)
    
    Returns:
        Protein sequence (single-letter amino acid codes)
    
    Example:
        >>> translate_sequence("ATGAAATAG")
        'MK*'
    """
    sequence = sequence.upper()
    protein = []
    
    for i in range(frame, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'X')
            protein.append(aa)
    
    return ''.join(protein)


def reverse_complement(sequence: str) -> str:
    """
    Generate reverse complement of DNA sequence.
    
    Args:
        sequence: DNA sequence string
    
    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    sequence = sequence.upper()
    return ''.join(complement.get(b, b) for b in reversed(sequence))


def calculate_gc(sequence: str) -> float:
    """Calculate GC content percentage."""
    sequence = sequence.upper()
    gc = sequence.count('G') + sequence.count('C')
    total = len([b for b in sequence if b in 'ATCG'])
    return (gc / total * 100) if total > 0 else 0.0


def find_orfs(
    sequence: str,
    min_length: int = 100,
    include_reverse: bool = True,
    use_alternative_starts: bool = False
) -> List[ORF]:
    """
    Find all Open Reading Frames in a DNA sequence.
    
    Scans all three reading frames on the forward strand,
    and optionally the reverse complement strand.
    
    Args:
        sequence: DNA sequence string
        min_length: Minimum ORF length in nucleotides (default 100)
        include_reverse: Whether to scan reverse complement (default True)
        use_alternative_starts: Include GTG, TTG, CTG start codons
    
    Returns:
        List of ORF objects sorted by length (descending)
    
    Example:
        >>> orfs = find_orfs("ATGAAAAAAAAATGA", min_length=9)
        >>> print(len(orfs), orfs[0].length_nt if orfs else 0)
    """
    sequence = sequence.upper()
    orfs = []
    
    start_codons = STANDARD_START_CODONS.copy()
    if use_alternative_starts:
        start_codons.update(ALTERNATIVE_START_CODONS)
    
    strands = [('+', sequence)]
    if include_reverse:
        strands.append(('-', reverse_complement(sequence)))
    
    for strand, seq in strands:
        for frame in range(3):
            orfs.extend(_find_orfs_in_frame(
                seq, frame, strand, start_codons, min_length
            ))
    
    orfs.sort(key=lambda x: x.length_nt, reverse=True)
    logger.info(f"Found {len(orfs)} ORFs >= {min_length} nt")
    return orfs


def _find_orfs_in_frame(
    sequence: str,
    frame: int,
    strand: str,
    start_codons: Set[str],
    min_length: int
) -> List[ORF]:
    """Find ORFs in a single reading frame."""
    orfs = []
    
    for i in range(frame, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        
        if codon in start_codons:
            for j in range(i + 3, len(sequence) - 2, 3):
                stop_codon = sequence[j:j+3]
                
                if stop_codon in STOP_CODONS:
                    orf_length = j + 3 - i
                    
                    if orf_length >= min_length:
                        orf_seq = sequence[i:j+3]
                        protein = translate_sequence(orf_seq)
                        
                        orfs.append(ORF(
                            start=i + 1,
                            end=j + 3,
                            length_nt=orf_length,
                            length_aa=len(protein) - 1,
                            frame=frame + 1,
                            strand=strand,
                            start_codon=codon,
                            stop_codon=stop_codon,
                            sequence=orf_seq,
                            protein=protein,
                            gc_content=calculate_gc(orf_seq)
                        ))
                    break
    
    return orfs


def get_orf_summary(orfs: List[ORF]) -> Dict[str, Any]:
    """
    Generate summary statistics for a list of ORFs.
    
    Args:
        orfs: List of ORF objects
    
    Returns:
        Dictionary with ORF statistics
    """
    if not orfs:
        return {'total': 0}
    
    lengths = [orf.length_nt for orf in orfs]
    frames = {}
    for orf in orfs:
        key = f"{orf.strand}{orf.frame}"
        frames[key] = frames.get(key, 0) + 1
    
    return {
        'total': len(orfs),
        'average_length': sum(lengths) / len(lengths),
        'max_length': max(lengths),
        'min_length': min(lengths),
        'by_frame': frames,
        'longest_orf': orfs[0] if orfs else None
    }
