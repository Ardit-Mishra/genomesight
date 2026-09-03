from collections import Counter, defaultdict
from typing import Dict

from Bio.Data import CodonTable

CODON_TO_AMINO_ACID = CodonTable.unambiguous_dna_by_name["Standard"].forward_table


def calculate_codon_usage(sequence: str) -> Dict[str, int]:
    """Calculate in-frame counts for standard sense codons."""
    sequence = sequence.upper().replace("U", "T")
    codons = (sequence[index:index + 3] for index in range(0, len(sequence) - 2, 3))
    return dict(Counter(codon for codon in codons if codon in CODON_TO_AMINO_ACID))


def calculate_rscu(sequence: str) -> Dict[str, float]:
    """Calculate RSCU within each amino acid's synonymous-codon family."""
    counts = calculate_codon_usage(sequence)
    families: dict[str, list[str]] = defaultdict(list)
    for codon, amino_acid in CODON_TO_AMINO_ACID.items():
        families[amino_acid].append(codon)

    rscu: Dict[str, float] = {}
    for codons in families.values():
        observed_total = sum(counts.get(codon, 0) for codon in codons)
        if not observed_total:
            continue
        expected_count = observed_total / len(codons)
        for codon in codons:
            if codon in counts:
                rscu[codon] = round(counts[codon] / expected_count, 3)
    return rscu
