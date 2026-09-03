import re

VALID_NUCLEOTIDES = frozenset("ACGTNRYSWKMBDHV")


def normalize_sequence(raw: str) -> str:
    """Normalize whitespace and RNA uracil without silently removing bases."""
    if not isinstance(raw, str) or not raw.strip():
        raise ValueError("A nucleotide sequence is required.")

    sequence = re.sub(r"\s+", "", raw).upper().replace("U", "T")
    invalid = sorted(set(sequence) - VALID_NUCLEOTIDES)
    if invalid:
        raise ValueError(
            "Unsupported nucleotide symbol(s): "
            f"{', '.join(invalid)}. Use DNA/RNA IUPAC nucleotide symbols only."
        )
    return sequence


def is_valid_dna(seq: str, allow_ambiguity: bool = True) -> bool:
    try:
        normalized = normalize_sequence(seq)
    except ValueError:
        return False
    return allow_ambiguity or set(normalized).issubset({"A", "C", "G", "T"})


def sanitize_raw_sequence(raw: str) -> str | None:
    try:
        return normalize_sequence(raw)
    except ValueError:
        return None


def wrap_raw_as_fasta(raw: str, header: str = ">sequence") -> str:
    seq = normalize_sequence(raw)
    lines = [header]
    for i in range(0, len(seq), 60):
        lines.append(seq[i:i + 60])
    return "\n".join(lines)
