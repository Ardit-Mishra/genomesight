"""
Codon Usage Analysis Module

This module counts codon usage in nucleotide sequences and computes Relative
Synonymous Codon Usage (RSCU) — the standard measure of whether an organism
(or a single ORF) prefers one synonymous codon over another for the same
amino acid.

Design Notes:
    - Reuses the standard genetic code from app.core.orf so the codon table
      is defined in exactly one place.
    - A codon is classified into exactly one of three buckets: a clean
      A/T/C/G triplet (counted), a triplet containing an IUPAC ambiguity
      code or other non-ACGT character (counted separately as ambiguous,
      never silently dropped), or a trailing partial codon of 1-2 leftover
      bases at the end of the reading frame (counted separately, never
      rounded away). Every base in the input is accounted for in one of
      these three tallies.
    - RSCU for a codon is left out of the result (not reported as 0.0) when
      its amino acid family was never observed at all — a fabricated zero
      would be indistinguishable from "this specific codon was avoided,"
      which is a real finding the true absence-of-data case is not.

RSCU Definition:
    For an amino acid encoded by n synonymous codons, with the amino acid
    observed `family_total` times total:
        expected_per_codon = family_total / n
        RSCU(codon) = observed_count(codon) / expected_per_codon
    RSCU = 1.0 means the codon is used exactly as often as expected under
    uniform usage across its synonymous family; > 1.0 means it is favored,
    < 1.0 means it is avoided.

Functions:
    count_codons: Tally codon usage (and RSCU) for one sequence in one frame
    codon_usage_for_orf: Codon usage for a single detected ORF
    merge_codon_usage: Combine several CodonUsageResult tallies into one
    get_codon_table_rows: Flatten a result into rows for display
"""

from typing import Any, Dict, List, Optional
from collections import Counter
from dataclasses import dataclass, field
import logging

from app.core.orf import CODON_TABLE

logger = logging.getLogger(__name__)

VALID_BASES = set('ATCG')

# Group the 64 standard codons by the amino acid (or stop, '*') they encode,
# so RSCU has a synonymous-family denominator to divide into.
_AA_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)


@dataclass
class CodonUsageResult:
    """Codon usage tally for one sequence, one ORF, or a merged aggregate."""
    codon_counts: Dict[str, int] = field(default_factory=dict)
    amino_acid_counts: Dict[str, int] = field(default_factory=dict)
    rscu: Dict[str, float] = field(default_factory=dict)
    total_codons: int = 0
    counted_codons: int = 0
    ambiguous_codons: int = 0
    partial_trailing_bases: int = 0
    source_length: int = 0


def _amino_acid_totals(codon_counts: Dict[str, int]) -> Dict[str, int]:
    """Sum codon counts up to their amino acid. Every key in codon_counts
    is one of the 64 standard triplets, so every one has an entry in
    CODON_TABLE — nothing here can fall through unclassified."""
    totals: Counter = Counter()
    for codon, count in codon_counts.items():
        aa = CODON_TABLE.get(codon)
        if aa is not None:
            totals[aa] += count
    return dict(totals)


def _calculate_rscu(codon_counts: Dict[str, int], amino_acid_counts: Dict[str, int]) -> Dict[str, float]:
    """RSCU per codon. A codon's whole amino-acid family is skipped (left
    out of the returned dict, not set to 0.0) when that family was never
    observed — there is no usage pattern to report, only a lack of data."""
    rscu: Dict[str, float] = {}
    for aa, codons in _AA_TO_CODONS.items():
        family_total = amino_acid_counts.get(aa, 0)
        if family_total == 0:
            continue
        expected_per_codon = family_total / len(codons)
        for codon in codons:
            observed = codon_counts.get(codon, 0)
            rscu[codon] = observed / expected_per_codon
    return rscu


def count_codons(sequence: str, frame: int = 0) -> CodonUsageResult:
    """
    Count codon usage across a nucleotide sequence read in one frame.

    Walks non-overlapping triplets starting at `frame`. Each triplet is
    classified into exactly one bucket:
        - all four bases A/T/C/G -> counted in codon_counts
        - contains an IUPAC ambiguity code (or anything else) -> counted in
          ambiguous_codons, never guessed at or dropped
    Any 1-2 bases left over after the last full triplet are counted in
    partial_trailing_bases rather than silently discarded.

    Args:
        sequence: Nucleotide sequence (DNA). Case-insensitive.
        frame: Reading frame offset — must be 0, 1, or 2.

    Returns:
        CodonUsageResult with raw counts, per-amino-acid totals, and RSCU.

    Raises:
        ValueError: If frame is not 0, 1, or 2.

    Example:
        >>> r = count_codons("ATGATGATC")
        >>> r.codon_counts
        {'ATG': 2, 'ATC': 1}
    """
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0, 1, or 2, got {frame}")

    sequence = sequence.upper()
    n = len(sequence)

    codon_counts: Counter = Counter()
    ambiguous_codons = 0

    i = frame
    while i + 3 <= n:
        codon = sequence[i:i + 3]
        if all(base in VALID_BASES for base in codon):
            codon_counts[codon] += 1
        else:
            ambiguous_codons += 1
        i += 3

    partial_trailing_bases = n - i

    amino_acid_counts = _amino_acid_totals(codon_counts)
    rscu = _calculate_rscu(codon_counts, amino_acid_counts)
    counted_codons = sum(codon_counts.values())

    logger.info(
        f"Counted {counted_codons} codons, {ambiguous_codons} ambiguous, "
        f"{partial_trailing_bases} trailing base(s) (frame {frame})"
    )

    return CodonUsageResult(
        codon_counts=dict(codon_counts),
        amino_acid_counts=amino_acid_counts,
        rscu=rscu,
        total_codons=counted_codons + ambiguous_codons,
        counted_codons=counted_codons,
        ambiguous_codons=ambiguous_codons,
        partial_trailing_bases=partial_trailing_bases,
        source_length=n,
    )


def codon_usage_for_orf(orf) -> CodonUsageResult:
    """
    Codon usage for a single detected ORF.

    An ORF's `.sequence` (from app.core.orf.find_orfs) already starts at its
    start codon, so it is read in frame 0 relative to itself — there is no
    separate "frame" concept for an individual ORF.

    Args:
        orf: An app.core.orf.ORF instance (or anything with a `.sequence`
             attribute holding the ORF's nucleotide sequence).

    Returns:
        CodonUsageResult for that ORF's sequence.
    """
    return count_codons(orf.sequence, frame=0)


def merge_codon_usage(results: List[CodonUsageResult]) -> CodonUsageResult:
    """
    Combine several codon-usage tallies (e.g. one per ORF) into one
    aggregate table, recomputing RSCU over the combined counts rather than
    averaging per-tally RSCU values (which would weight a short ORF equally
    with a long one).

    Args:
        results: List of CodonUsageResult, e.g. from codon_usage_for_orf
                 applied to each detected ORF.

    Returns:
        A single merged CodonUsageResult. Returns an all-zero result (not
        an error) for an empty list — "no ORFs" is a valid, reportable
        state, not a failure.
    """
    codon_counts: Counter = Counter()
    ambiguous_codons = 0
    partial_trailing_bases = 0
    total_codons = 0
    source_length = 0

    for result in results:
        codon_counts.update(result.codon_counts)
        ambiguous_codons += result.ambiguous_codons
        partial_trailing_bases += result.partial_trailing_bases
        total_codons += result.total_codons
        source_length += result.source_length

    amino_acid_counts = _amino_acid_totals(codon_counts)
    rscu = _calculate_rscu(codon_counts, amino_acid_counts)
    counted_codons = sum(codon_counts.values())

    return CodonUsageResult(
        codon_counts=dict(codon_counts),
        amino_acid_counts=amino_acid_counts,
        rscu=rscu,
        total_codons=total_codons,
        counted_codons=counted_codons,
        ambiguous_codons=ambiguous_codons,
        partial_trailing_bases=partial_trailing_bases,
        source_length=source_length,
    )


def get_codon_table_rows(result: CodonUsageResult) -> List[Dict[str, Any]]:
    """
    Flatten a CodonUsageResult into rows for a UI table: one row per one of
    the 64 standard codons, sorted by amino acid then by descending count.

    `fraction_of_aa_pct` and `rscu` are None (not 0) for a codon whose whole
    amino-acid family was never observed — there is nothing to report a
    fraction or preference of.

    Args:
        result: A CodonUsageResult from count_codons/codon_usage_for_orf/
                merge_codon_usage.

    Returns:
        List of dicts with keys: codon, amino_acid, count,
        fraction_of_aa_pct, rscu.
    """
    rows = []
    for codon in sorted(CODON_TABLE.keys()):
        aa = CODON_TABLE[codon]
        count = result.codon_counts.get(codon, 0)
        family_total = result.amino_acid_counts.get(aa, 0)
        fraction = (count / family_total * 100) if family_total > 0 else None
        rows.append({
            'codon': codon,
            'amino_acid': aa,
            'count': count,
            'fraction_of_aa_pct': fraction,
            'rscu': result.rscu.get(codon),
        })

    rows.sort(key=lambda r: (r['amino_acid'], -r['count']))
    return rows
