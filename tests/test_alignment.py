"""
Tests for pairwise sequence alignment.

The module existed and was never reachable from the UI, so nothing exercised it
and a parsing bug went unnoticed: `_format_alignment` read Biopython's
`alignment.format()` output, which since Biopython 1.80 includes coordinate
columns ("target  0 ATGC... 14"). Those labels were consumed as sequence, so two
14 bp sequences differing at one position reported 81.1% identity over a
37-position alignment instead of 92.9% over 14.

These tests pin the arithmetic to the sequences themselves.
"""

import pytest

from app.core.alignment import (
    AlignmentResult,
    align_sequences,
    calculate_identity,
    format_alignment_display,
)


class TestAlignmentArithmetic:
    """Identity, length and gap counts must describe the sequences, not the display."""

    def test_single_mismatch(self):
        result = align_sequences("ATGCGTACGTTAGC", "ATGCGTTCGTTAGC")

        # 14 bp each, differing at one position, no indels needed.
        assert result.alignment_length == 14
        assert result.identity_count == 13
        assert result.gaps == 0
        assert result.identity == pytest.approx(92.857, abs=0.01)

    def test_identical_sequences_are_100_percent(self):
        seq = "ATGCGTACGTTAGC"
        result = align_sequences(seq, seq)

        assert result.identity == pytest.approx(100.0)
        assert result.identity_count == result.alignment_length == len(seq)
        assert result.gaps == 0

    def test_deletion_introduces_gaps(self):
        result = align_sequences("ATGCGTACGTTAGC", "ATGCGTTAGC")

        assert result.gaps > 0
        assert result.alignment_length >= 14
        assert 0 < result.identity < 100

    def test_alignment_rows_contain_only_sequence_characters(self):
        """The regression guard: no labels, digits or coordinates in the rows."""
        result = align_sequences("ATGCGTACGTTAGC", "ATGCGTTCGTTAGC")

        for row in (result.aligned_seq1, result.aligned_seq2):
            assert row, "aligned row should not be empty"
            assert set(row) <= set("ACGTUN-"), f"unexpected characters in {row!r}"
            assert "target" not in row and "query" not in row

    def test_rows_are_equal_length(self):
        """Aligned rows are padded to the same length; identity depends on it."""
        result = align_sequences("ATGCGTACGTTAGC", "ATGCGTTAGC")
        assert len(result.aligned_seq1) == len(result.aligned_seq2)
        assert len(result.match_line) == len(result.aligned_seq1)


class TestAlignmentModes:
    def test_global_and_local_both_run(self):
        for mode in ("global", "local"):
            result = align_sequences("ATGCGTACGTTAGC", "CGTACG", mode=mode)
            assert isinstance(result, AlignmentResult)
            assert result.mode == mode
            assert result.alignment_length > 0

    def test_local_finds_the_embedded_subsequence(self):
        """A local alignment of an exact substring should match it perfectly."""
        result = align_sequences("TTTTATGCGTACGTTAGCTTTT", "ATGCGTACGTTAGC", mode="local")
        assert result.identity == pytest.approx(100.0)


class TestEdgeCases:
    def test_empty_sequence_returns_empty_result(self):
        result = align_sequences("", "ATGC")
        assert result.alignment_length == 0
        assert result.identity == 0

    def test_lowercase_and_spaces_are_normalised(self):
        spaced = align_sequences("atg cgt acg tta gc", "ATGCGTACGTTAGC")
        assert spaced.identity == pytest.approx(100.0)


class TestDisplay:
    def test_display_reports_the_same_numbers_as_the_result(self):
        result = align_sequences("ATGCGTACGTTAGC", "ATGCGTTCGTTAGC")
        text = format_alignment_display(result)

        assert "92.9%" in text
        assert "(13/14)" in text
        assert "Global" in text

    def test_calculate_identity_matches_alignment_identity(self):
        a, b = "ATGCGTACGTTAGC", "ATGCGTTCGTTAGC"
        assert calculate_identity(a, b) == pytest.approx(align_sequences(a, b).identity, abs=0.01)


class TestFailureIsDistinguishableFromAResult:
    """A failed alignment and a genuine zero-similarity alignment used to be
    byte-identical: 0.0% identity, score 0.0, 0 gaps, empty strings. The UI
    rendered both as four confident metrics, so a crashed aligner was
    indistinguishable from a real finding about two unrelated sequences.
    """

    def test_successful_alignment_is_marked_ok_and_carries_no_error(self):
        result = align_sequences("ATGCGTACGTTAGC", "ATGCGTTCGTTAGC")
        assert result.ok is True
        assert result.error is None

    def test_failure_is_flagged_rather_than_returned_as_zero_identity(self):
        # An unsupported mode raises inside Biopython. Before the fix this
        # escaped the function entirely and took the page down with it.
        result = align_sequences("ATGC", "ATGC", mode="not-a-real-mode")

        assert result.ok is False
        assert result.error is not None
        assert "invalid mode" in result.error.lower()

    def test_a_real_zero_similarity_result_is_not_flagged_as_an_error(self):
        # The distinction only matters if genuine dissimilarity stays clean.
        result = align_sequences("AAAAAAAA", "TTTTTTTT", mode="local")
        assert result.ok is True
        assert result.error is None
