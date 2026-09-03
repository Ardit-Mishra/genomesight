"""
Motif Discovery and Pattern Matching Module

This module provides functions for searching DNA sequences for specific patterns
using IUPAC ambiguity codes. It includes a library of common restriction enzyme
recognition sites.

IUPAC Ambiguity Codes:
    R = A or G (purine)
    Y = C or T (pyrimidine)
    S = G or C (strong)
    W = A or T (weak)
    K = G or T (keto)
    M = A or C (amino)
    B = C, G, or T (not A)
    D = A, G, or T (not C)
    H = A, C, or T (not G)
    V = A, C, or G (not T)
    N = A, C, G, or T (any)

Design Notes:
    - Converts IUPAC patterns to regex for flexible matching
    - Reports both 0-based and 1-based positions
    - Includes highlighted sequence output for visualization
"""

import re
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


IUPAC_CODES = {
    'A': 'A',
    'T': 'T',
    'C': 'C',
    'G': 'G',
    'U': 'U',
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGT]',
    'D': '[AGT]',
    'H': '[ACT]',
    'V': '[ACG]',
    'N': '[ACGT]'
}


RESTRICTION_ENZYMES = {
    'EcoRI': 'GAATTC',
    'BamHI': 'GGATCC',
    'HindIII': 'AAGCTT',
    'EcoRV': 'GATATC',
    'PstI': 'CTGCAG',
    'SalI': 'GTCGAC',
    'XbaI': 'TCTAGA',
    'NotI': 'GCGGCCGC',
    'SmaI': 'CCCGGG',
    'KpnI': 'GGTACC'
}


@dataclass
class MotifMatch:
    """Data class representing a motif match."""
    pattern: str
    sequence_id: str
    start_0based: int
    start_1based: int
    end: int
    matched_sequence: str
    context: str


def iupac_to_regex(pattern: str) -> str:
    """
    Convert an IUPAC pattern to a regular expression.
    
    Args:
        pattern: IUPAC nucleotide pattern (e.g., "ATGR" for ATG[AG])
    
    Returns:
        Regular expression string
    
    Example:
        >>> iupac_to_regex("ATGR")
        'ATG[AG]'
    """
    pattern = pattern.upper()
    regex_parts = []
    
    for char in pattern:
        if char in IUPAC_CODES:
            regex_parts.append(IUPAC_CODES[char])
        else:
            regex_parts.append(re.escape(char))
    
    return ''.join(regex_parts)


def search_motif(
    sequences: List,
    pattern: str,
    context_length: int = 10
) -> List[MotifMatch]:
    """
    Search for a motif pattern in sequences using IUPAC codes.
    
    Args:
        sequences: List of Bio.SeqRecord objects
        pattern: IUPAC nucleotide pattern to search for
        context_length: Bases to include before/after match for context
    
    Returns:
        List of MotifMatch objects
    
    Example:
        >>> matches = search_motif(sequences, "ATGR")
        >>> for m in matches:
        ...     print(f"Found at position {m.start_1based}")
    """
    if not pattern:
        return []
    
    regex_pattern = iupac_to_regex(pattern.upper())
    compiled_pattern = re.compile(regex_pattern)
    matches = []
    
    for i, seq in enumerate(sequences):
        sequence_str = str(seq.seq).upper()
        seq_id = seq.id if hasattr(seq, 'id') else f"Sequence_{i}"
        
        for match in compiled_pattern.finditer(sequence_str):
            start = match.start()
            end = match.end()
            
            context_start = max(0, start - context_length)
            context_end = min(len(sequence_str), end + context_length)
            context = sequence_str[context_start:context_end]
            
            matches.append(MotifMatch(
                pattern=pattern.upper(),
                sequence_id=seq_id,
                start_0based=start,
                start_1based=start + 1,
                end=end,
                matched_sequence=match.group(),
                context=context
            ))
    
    logger.info(f"Found {len(matches)} matches for pattern '{pattern}'")
    return matches


def search_restriction_site(
    sequences: List,
    enzyme: str
) -> List[MotifMatch]:
    """
    Search for a restriction enzyme recognition site.
    
    Args:
        sequences: List of Bio.SeqRecord objects
        enzyme: Name of restriction enzyme (e.g., "EcoRI")
    
    Returns:
        List of MotifMatch objects
    
    Raises:
        ValueError: If enzyme not in database
    """
    if enzyme not in RESTRICTION_ENZYMES:
        available = ', '.join(sorted(RESTRICTION_ENZYMES.keys()))
        raise ValueError(f"Unknown enzyme: {enzyme}. Available: {available}")
    
    pattern = RESTRICTION_ENZYMES[enzyme]
    return search_motif(sequences, pattern)


def highlight_matches(
    sequence: str,
    matches: List[MotifMatch],
    highlight_start: str = "[",
    highlight_end: str = "]"
) -> str:
    """
    Create a highlighted version of sequence showing matches.
    
    Args:
        sequence: Original sequence string
        matches: List of MotifMatch objects for this sequence
        highlight_start: Character(s) to insert before match
        highlight_end: Character(s) to insert after match
    
    Returns:
        Sequence with matches highlighted
    
    Example:
        >>> highlight_matches("ATGCATGC", matches)
        "AT[GC]ATGC"  # if GC was matched
    """
    if not matches:
        return sequence
    
    sorted_matches = sorted(matches, key=lambda m: m.start_0based, reverse=True)
    
    result = sequence
    for match in sorted_matches:
        start = match.start_0based
        end = match.end
        result = (
            result[:start] + 
            highlight_start + 
            result[start:end] + 
            highlight_end + 
            result[end:]
        )
    
    return result


def find_common_motifs(
    sequences: List,
    min_length: int = 6,
    max_length: int = 12,
    min_conservation: float = 0.3
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Discover conserved motifs across multiple sequences.
    
    Uses k-mer analysis to find patterns present in multiple sequences.
    
    Args:
        sequences: List of Bio.SeqRecord objects
        min_length: Minimum motif length
        max_length: Maximum motif length
        min_conservation: Minimum fraction of sequences containing motif
    
    Returns:
        Dictionary mapping k-mer size to list of enriched motifs
    """
    from collections import Counter
    
    results = {}
    
    for k in range(min_length, max_length + 1):
        kmer_counts = Counter()
        total_kmers = 0
        
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i+k]
                if all(base in 'ATCG' for base in kmer):
                    kmer_counts[kmer] += 1
                    total_kmers += 1
        
        enriched = []
        for kmer, count in kmer_counts.most_common(20):
            seq_count = sum(
                1 for seq in sequences 
                if kmer in str(seq.seq).upper()
            )
            conservation = seq_count / len(sequences) if sequences else 0
            
            if conservation >= min_conservation:
                enriched.append({
                    'motif': kmer,
                    'count': count,
                    'frequency': count / total_kmers if total_kmers > 0 else 0,
                    'conservation': conservation,
                    'sequences_with_motif': seq_count
                })
        
        if enriched:
            results[f'{k}-mer'] = enriched[:10]
    
    return results


def get_enzyme_list() -> Dict[str, str]:
    """Return dictionary of available restriction enzymes."""
    return RESTRICTION_ENZYMES.copy()
