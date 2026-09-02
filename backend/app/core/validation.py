from Bio.Seq import Seq
from typing import Dict, List, Optional

VALID_NUCLEOTIDES = set("ACGTNRYSWKMBDHV")

def is_valid_dna(seq: str, allow_ambiguity: bool = True) -> bool:
    if not seq or not isinstance(seq, str):
        return False
    allowed = VALID_NUCLEOTIDES if allow_ambiguity else set("ACGT")
    return all(base.upper() in allowed for base in seq)

def sanitize_raw_sequence(raw: str) -> Optional[str]:
    if not raw or not isinstance(raw, str):
        return None
    cleaned = ''.join(ch.upper() for ch in raw if ch.isalpha())
    # Allow only standard DNA + common ambiguity
    allowed_chars = set("ACGTNRYSWKMBDHV")
    cleaned = ''.join(ch for ch in cleaned if ch in allowed_chars)
    if len(cleaned) < 3:
        return None
    return cleaned

def wrap_raw_as_fasta(raw: str, header: str = ">sequence") -> str:
    seq = sanitize_raw_sequence(raw)
    if seq is None:
        return ""
    # Wrap at 60 chars for readability
    lines = [header]
    for i in range(0, len(seq), 60):
        lines.append(seq[i:i+60])
    return "\n".join(lines)
