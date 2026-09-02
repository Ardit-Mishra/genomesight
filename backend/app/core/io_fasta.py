from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import io
from typing import List

def parse_fasta_string(content: str) -> List[SeqRecord]:
    """Parse FASTA or FASTQ content string into a list of SeqRecord objects."""
    if not content.strip():
        return []

    file_like = io.StringIO(content.strip())
    format_type = "fastq" if content.strip().startswith("@") else "fasta"

    try:
        records = list(SeqIO.parse(file_like, format_type))
        return records
    except Exception:
        file_like.seek(0)
        return list(SeqIO.parse(file_like, "fasta"))
