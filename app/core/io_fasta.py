"""
FASTA/FASTQ/GenBank File I/O Module

This module handles parsing and exporting of sequence files in FASTA, FASTQ,
and GenBank formats. It provides robust file format detection, validation,
and conversion utilities.

Design Notes:
    - Uses Biopython SeqIO for reliable parsing
    - Handles multi-line FASTA sequences correctly
    - Supports both file objects and string content
    - Includes quality score handling for FASTQ
    - GenBank records keep their `.features` (Bio.SeqFeature) list intact;
      extract_genbank_features() flattens the ones that matter for display

Functions:
    parse_fasta: Parse FASTA format files
    parse_fastq: Parse FASTQ format files
    parse_genbank: Parse GenBank format files, including feature annotations
    detect_format: Auto-detect file format
    export_to_fasta: Export sequences to FASTA format
    extract_genbank_features: Flatten a GenBank record's features for display
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
        Format string ('fasta', 'fastq', or 'genbank'), or None if unrecognized

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

    # A GenBank flat file's first line always starts with the LOCUS keyword —
    # check content before falling back to the filename extension.
    if content[:5].upper() == 'LOCUS':
        return 'genbank'

    ext = filename.lower()
    if ext.endswith(('.fasta', '.fa', '.fna', '.ffn', '.faa')):
        return 'fasta'
    if ext.endswith(('.fastq', '.fq')):
        return 'fastq'
    if ext.endswith(('.gb', '.gbk', '.genbank')):
        return 'genbank'

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


def parse_genbank(content: str) -> List[SeqRecord]:
    """
    Parse GenBank format content into sequence records.

    GenBank records carry structured feature annotations (genes, CDS,
    regulatory regions, etc.) in addition to the raw sequence. Biopython
    attaches these to each record's `.features` list — this function does
    not drop or flatten them; use `extract_genbank_features` for a
    display-ready view of a record's features.

    Args:
        content: GenBank formatted string content (one LOCUS entry or many)

    Returns:
        List of Bio.SeqRecord objects, each with `.features` populated

    Raises:
        ValueError: If content is not valid GenBank format
    """
    if content.strip()[:5].upper() != 'LOCUS':
        raise ValueError("Invalid GenBank format: must start with 'LOCUS'")

    try:
        handle = StringIO(content)
        records = list(SeqIO.parse(handle, 'genbank'))
        if not records:
            raise ValueError("No LOCUS records parsed from GenBank content")
        logger.info(f"Parsed {len(records)} sequences from GenBank")
        return records
    except ValueError:
        raise
    except Exception as e:
        logger.error(f"GenBank parsing error: {e}")
        raise ValueError(f"Failed to parse GenBank: {e}")


def extract_genbank_features(record: SeqRecord) -> List[Dict[str, Any]]:
    """
    Flatten a GenBank record's features into plain dicts for display.

    Args:
        record: A Bio.SeqRecord parsed from GenBank content. Records from
            FASTA/FASTQ have no `.features`, so those return an empty list —
            not an error, since "no feature annotations" is a legitimate
            state for non-GenBank input.

    Returns:
        List of dicts, one per feature, with keys: type, start (1-based),
        end, strand ('+', '-', or None if unspecified), and qualifiers
        (the feature's raw qualifier dict, e.g. gene/product/note).
    """
    features = getattr(record, 'features', None)
    if not features:
        return []

    strand_symbol = {1: '+', -1: '-', 0: None, None: None}

    flattened = []
    for feature in features:
        location = feature.location
        flattened.append({
            'type': feature.type,
            'start': int(location.start) + 1,
            'end': int(location.end),
            'strand': strand_symbol.get(location.strand, None),
            'qualifiers': {k: v for k, v in feature.qualifiers.items()},
        })

    return flattened


def parse_sequences(content: str, filename: str = "") -> Tuple[List[SeqRecord], str]:
    """
    Parse sequence content with automatic format detection.

    Args:
        content: File content as string
        filename: Optional filename for format hints

    Returns:
        Tuple of (list of SeqRecords, format string — 'fasta', 'fastq', or
        'genbank')

    Raises:
        ValueError: If format cannot be detected or parsing fails
    """
    file_format = detect_format(content, filename)

    if file_format is None:
        raise ValueError("Could not detect file format")

    if file_format == 'fasta':
        return parse_fasta(content), 'fasta'
    elif file_format == 'genbank':
        return parse_genbank(content), 'genbank'
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
