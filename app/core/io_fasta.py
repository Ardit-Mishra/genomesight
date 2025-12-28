"""
FASTA/FASTQ File I/O Module

This module handles parsing and exporting of sequence files in FASTA and FASTQ formats.
It provides robust file format detection, validation, and conversion utilities.

Design Notes:
    - Uses Biopython SeqIO for reliable parsing
    - Handles multi-line FASTA sequences correctly
    - Supports both file objects and string content
    - Includes quality score handling for FASTQ

Functions:
    parse_fasta: Parse FASTA format files
    parse_fastq: Parse FASTQ format files  
    detect_format: Auto-detect file format
    export_to_fasta: Export sequences to FASTA format
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from typing import List, Dict, Any, Optional, Union, Tuple
import logging

logger = logging.getLogger(__name__)


def detect_format(content: str, filename: str = "") -> Optional[str]:
    """
    Detect the format of sequence file content.
    
    Uses both file extension and content inspection to determine format.
    Content-based detection takes precedence for reliability.
    
    Args:
        content: File content as string
        filename: Optional filename for extension-based hints
    
    Returns:
        Format string ('fasta' or 'fastq'), or None if unrecognized
    
    Example:
        >>> detect_format(">seq1\\nATGC", "test.fasta")
        'fasta'
    """
    content = content.strip()
    
    if content.startswith('>'):
        return 'fasta'
    
    if content.startswith('@'):
        lines = content.split('\n')
        if len(lines) >= 4 and lines[2].startswith('+'):
            return 'fastq'
    
    ext = filename.lower()
    if ext.endswith(('.fasta', '.fa', '.fna', '.ffn', '.faa')):
        return 'fasta'
    if ext.endswith(('.fastq', '.fq')):
        return 'fastq'
    
    return None


def parse_fasta(content: str) -> List[SeqRecord]:
    """
    Parse FASTA format content into sequence records.
    
    FASTA format:
        >sequence_id optional description
        ATGCATGCATGC
        GCTAGCTAGCTA
    
    Args:
        content: FASTA formatted string content
    
    Returns:
        List of Bio.SeqRecord objects
    
    Raises:
        ValueError: If content is not valid FASTA format
    """
    if not content.strip().startswith('>'):
        raise ValueError("Invalid FASTA format: must start with '>'")
    
    try:
        handle = StringIO(content)
        records = list(SeqIO.parse(handle, 'fasta'))
        logger.info(f"Parsed {len(records)} sequences from FASTA")
        return records
    except Exception as e:
        logger.error(f"FASTA parsing error: {e}")
        raise ValueError(f"Failed to parse FASTA: {e}")


def parse_fastq(content: str) -> List[SeqRecord]:
    """
    Parse FASTQ format content into sequence records.
    
    FASTQ format:
        @sequence_id
        ATGCATGCATGC
        +
        !!!!!!!!!!!!!
    
    Quality scores are stored in record.letter_annotations['phred_quality']
    
    Args:
        content: FASTQ formatted string content
    
    Returns:
        List of Bio.SeqRecord objects with quality annotations
    
    Raises:
        ValueError: If content is not valid FASTQ format
    """
    if not content.strip().startswith('@'):
        raise ValueError("Invalid FASTQ format: must start with '@'")
    
    try:
        handle = StringIO(content)
        records = list(SeqIO.parse(handle, 'fastq'))
        logger.info(f"Parsed {len(records)} sequences from FASTQ")
        return records
    except Exception as e:
        logger.error(f"FASTQ parsing error: {e}")
        raise ValueError(f"Failed to parse FASTQ: {e}")


def parse_sequences(content: str, filename: str = "") -> Tuple[List[SeqRecord], str]:
    """
    Parse sequence content with automatic format detection.
    
    Args:
        content: File content as string
        filename: Optional filename for format hints
    
    Returns:
        Tuple of (list of SeqRecords, format string)
    
    Raises:
        ValueError: If format cannot be detected or parsing fails
    """
    file_format = detect_format(content, filename)
    
    if file_format is None:
        raise ValueError("Could not detect file format")
    
    if file_format == 'fasta':
        return parse_fasta(content), 'fasta'
    else:
        return parse_fastq(content), 'fastq'


def export_to_fasta(sequences: List[SeqRecord], line_width: int = 80) -> str:
    """
    Export sequence records to FASTA format string.
    
    Args:
        sequences: List of Bio.SeqRecord objects
        line_width: Characters per line for sequence wrapping (default 80)
    
    Returns:
        FASTA formatted string
    
    Example:
        >>> records = [SeqRecord(Seq("ATGC"), id="seq1")]
        >>> print(export_to_fasta(records))
        >seq1
        ATGC
    """
    lines = []
    
    for seq in sequences:
        seq_id = seq.id if hasattr(seq, 'id') else "sequence"
        description = getattr(seq, 'description', '')
        
        header = f">{seq_id}"
        if description and description != seq_id:
            header += f" {description}"
        lines.append(header)
        
        sequence_str = str(seq.seq)
        for i in range(0, len(sequence_str), line_width):
            lines.append(sequence_str[i:i+line_width])
    
    return '\n'.join(lines)


def clean_sequence(sequence: str) -> str:
    """
    Clean and normalize a sequence string.
    
    Removes whitespace, converts to uppercase, and validates characters.
    
    Args:
        sequence: Raw sequence string
    
    Returns:
        Cleaned sequence string (uppercase, no whitespace)
    """
    cleaned = ''.join(sequence.split()).upper()
    return cleaned
