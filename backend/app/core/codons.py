from collections import Counter
from typing import Dict

def calculate_codon_usage(sequence: str) -> Dict[str, int]:
    """Calculate codon counts in a nucleotide sequence."""
    seq = sequence.upper().replace("U", "T")
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    return dict(Counter(codons))

def calculate_rscu(sequence: str) -> Dict[str, float]:
    """Calculate Relative Synonymous Codon Usage (RSCU)."""
    counts = calculate_codon_usage(sequence)
    total = sum(counts.values())
    if total == 0:
        return {}
    return {k: v / total for k, v in counts.items()}
