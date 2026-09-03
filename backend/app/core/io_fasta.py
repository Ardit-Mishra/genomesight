"""Sequence file parsing: FASTA, FASTQ and GenBank.

Format is decided by CONTENT, not by a filename that may be absent or wrong.
GenBank support is why this is not a two-line wrapper: a GenBank flat file
starts with the LOCUS keyword and parses under a different SeqIO format, so a
"starts with > or @" test silently rejects it -- which is exactly what the
previous version did, while the project documented GenBank as supported.
"""
from typing import List, Optional

import io

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def detect_format(content: str, filename: str = "") -> Optional[str]:
    """Return 'fasta', 'fastq', 'genbank', or None.

    Content wins over the extension, because content is what actually gets parsed.
    """
    content = content.strip()
    if not content:
        return None

    if content.startswith(">"):
        return "fasta"

    if content.startswith("@"):
        lines = content.split("\n")
        # A FASTQ record is exactly four lines with '+' third. Without that check an
        # '@'-prefixed FASTA-like file is mis-detected and parses to nothing.
        if len(lines) >= 4 and lines[2].startswith("+"):
            return "fastq"

    if content[:5].upper() == "LOCUS":
        return "genbank"

    ext = filename.lower()
    if ext.endswith((".fasta", ".fa", ".fna", ".ffn", ".faa")):
        return "fasta"
    if ext.endswith((".fastq", ".fq")):
        return "fastq"
    if ext.endswith((".gb", ".gbk", ".genbank")):
        return "genbank"

    return None


def parse_fasta_string(content: str, filename: str = "") -> List[SeqRecord]:
    """Parse FASTA, FASTQ or GenBank content into SeqRecords.

    Returns an empty list when nothing parses. The caller decides what an empty
    result means; this function does not invent records to avoid one.
    """
    if not content.strip():
        return []

    fmt = detect_format(content, filename)
    if fmt is None:
        # Unrecognised: try FASTA rather than refuse outright, since a bare
        # header-less sequence is wrapped upstream and still arrives as FASTA.
        fmt = "fasta"

    try:
        return list(SeqIO.parse(io.StringIO(content.strip()), fmt))
    except Exception:
        # A malformed file of a known format should not be reported as a
        # different format's failure; fall back once, then give up honestly.
        if fmt == "fasta":
            return []
        try:
            return list(SeqIO.parse(io.StringIO(content.strip()), "fasta"))
        except Exception:
            return []
