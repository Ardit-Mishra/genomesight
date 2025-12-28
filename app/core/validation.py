"""
Sequence Validation Module

This module provides functions for validating DNA/RNA sequences and
detecting sequence types. It includes IUPAC ambiguity code support.

Design Notes:
    - Supports standard IUPAC nucleotide codes
    - Distinguishes between DNA (contains T) and RNA (contains U)
    - Provides detailed validation feedback

IUPAC Nucleotide Codes:
    A, T, C, G - Standard bases
    U - Uracil (RNA)
    R - A or G (purine)
    Y - C or T (pyrimidine)
    S - G or C
    W - A or T
    K - G or T
    M - A or C
    B - C, G, or T
    D - A, G, or T
    H - A, C, or T
    V - A, C, or G
    N - Any base
"""

from typing import Dict, Any, List, Set, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


VALID_DNA_BASES = set('ATCGN')
VALID_RNA_BASES = set('AUCGN')
IUPAC_CODES = set('ATCGURYSWKMBDHVN')


@dataclass
class ValidationResult:
    """Result of sequence validation."""
    is_valid: bool
    sequence_type: str
    warnings: List[str]
    invalid_chars: Set[str]
    base_counts: Dict[str, int]


def validate_sequence(sequence: str) -> ValidationResult:
    """
    Validate a nucleotide sequence and detect its type.
    
    Checks for valid characters, detects DNA vs RNA,
    and provides warnings for ambiguous bases.
    
    Args:
        sequence: Nucleotide sequence string
    
    Returns:
        ValidationResult with validation details
    
    Example:
        >>> result = validate_sequence("ATGCATGC")
        >>> print(result.is_valid, result.sequence_type)
        True DNA
    """
    sequence = sequence.upper()
    
    base_counts = {}
    invalid_chars = set()
    warnings = []
    
    for char in sequence:
        if char.isspace():
            continue
        base_counts[char] = base_counts.get(char, 0) + 1
        if char not in IUPAC_CODES:
            invalid_chars.add(char)
    
    is_valid = len(invalid_chars) == 0
    
    seq_type = detect_sequence_type(sequence)
    
    ambiguous = set(sequence.upper()) & set('RYSWKMBDHVN')
    if ambiguous:
        warnings.append(f"Contains ambiguous bases: {', '.join(sorted(ambiguous))}")
    
    n_count = base_counts.get('N', 0)
    if n_count > 0:
        pct = (n_count / len(sequence)) * 100
        if pct > 5:
            warnings.append(f"High N content: {pct:.1f}%")
    
    if not is_valid:
        warnings.append(f"Invalid characters found: {', '.join(sorted(invalid_chars))}")
    
    return ValidationResult(
        is_valid=is_valid,
        sequence_type=seq_type,
        warnings=warnings,
        invalid_chars=invalid_chars,
        base_counts=base_counts
    )


def detect_sequence_type(sequence: str) -> str:
    """
    Determine if a sequence is DNA or RNA.
    
    DNA contains thymine (T), RNA contains uracil (U).
    If both are present, returns 'Mixed'.
    
    Args:
        sequence: Nucleotide sequence string
    
    Returns:
        'DNA', 'RNA', 'Mixed', or 'Unknown'
    """
    sequence = sequence.upper()
    has_t = 'T' in sequence
    has_u = 'U' in sequence
    
    if has_u and not has_t:
        return 'RNA'
    elif has_t and not has_u:
        return 'DNA'
    elif has_t and has_u:
        return 'Mixed'
    else:
        return 'DNA'


def validate_sequences(sequences: List) -> Dict[str, Any]:
    """
    Validate a list of sequence records.
    
    Provides aggregate validation results and recommendations.
    
    Args:
        sequences: List of Bio.SeqRecord objects
    
    Returns:
        Dictionary with validation summary and recommendations
    """
    if not sequences:
        return {
            'is_valid': False,
            'warnings': ['No sequences provided'],
            'recommendations': [],
            'quality_score': 0
        }
    
    results = {
        'is_valid': True,
        'warnings': [],
        'recommendations': [],
        'quality_score': 100
    }
    
    lengths = [len(seq.seq) for seq in sequences]
    avg_length = sum(lengths) / len(lengths)
    
    if avg_length < 10:
        results['warnings'].append("Very short sequences (avg < 10 bp)")
        results['recommendations'].append("Check if sequences are complete")
        results['quality_score'] -= 20
    
    all_invalid = set()
    for seq in sequences[:10]:
        validation = validate_sequence(str(seq.seq))
        all_invalid.update(validation.invalid_chars)
        if not validation.is_valid:
            results['is_valid'] = False
    
    if all_invalid:
        results['warnings'].append(f"Invalid characters: {', '.join(sorted(all_invalid))}")
        results['recommendations'].append("Consider cleaning sequence data")
        results['quality_score'] -= 10
    
    if len(sequences) < 5:
        results['recommendations'].append("Small dataset - consider adding more sequences")
    
    if len(sequences) > 10000:
        results['recommendations'].append("Large dataset - analysis may take longer")
    
    return results


def get_sequence_summary(sequence: str) -> Dict[str, Any]:
    """
    Generate a summary of sequence properties.
    
    Args:
        sequence: Nucleotide sequence string
    
    Returns:
        Dictionary with length, type, and base composition
    """
    validation = validate_sequence(sequence)
    
    total = sum(validation.base_counts.values())
    composition_pct = {
        base: (count / total) * 100 
        for base, count in validation.base_counts.items()
    } if total > 0 else {}
    
    return {
        'length': len(sequence.replace(' ', '')),
        'type': validation.sequence_type,
        'is_valid': validation.is_valid,
        'base_counts': validation.base_counts,
        'composition_percent': composition_pct,
        'warnings': validation.warnings
    }
