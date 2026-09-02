from collections import Counter
from typing import List, Tuple, Dict

def count_kmers(sequences: List[str], k: int = 3) -> Tuple[Dict[str, int], str]:
    """Count k-mers across sequences with pure Python fallback."""
    counts = Counter()
    for seq in sequences:
        seq = seq.upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if len(kmer) == k:
                counts[kmer] += 1
    return dict(counts), "python"
