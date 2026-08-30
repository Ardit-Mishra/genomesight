"""
Unit tests for the codon usage module.

Run with: pytest tests/test_codons.py -v
"""

import pytest

from app.core.codons import (
    count_codons,
    codon_usage_for_orf,
    merge_codon_usage,
    get_codon_table_rows,
    CodonUsageResult,
)
from app.core.orf import find_orfs


class TestBasicCounting:
    """Counts must land in exactly one bucket, and nothing should vanish."""

    def test_simple_counts(self):
        result = count_codons("ATGATGATC")
        assert result.codon_counts == {"ATG": 2, "ATC": 1}
        assert result.counted_codons == 3
        assert result.ambiguous_codons == 0
        assert result.partial_trailing_bases == 0
        assert result.total_codons == 3

    def test_frame_offset_changes_the_triplets(self):
        # "AATGATGATC" read from frame 0 vs frame 1 tiles completely
        # differently, so the two tallies must differ.
        seq = "AATGATGATC"
        frame0 = count_codons(seq, frame=0)
        frame1 = count_codons(seq, frame=1)
        assert frame0.codon_counts != frame1.codon_counts

    def test_invalid_frame_raises(self):
        with pytest.raises(ValueError):
            count_codons("ATGATG", frame=3)

    def test_lowercase_is_normalised(self):
        upper = count_codons("ATGATC")
        lower = count_codons("atgatc")
        assert upper.codon_counts == lower.codon_counts

    def test_empty_sequence(self):
        result = count_codons("")
        assert result.codon_counts == {}
        assert result.total_codons == 0
        assert result.partial_trailing_bases == 0


class TestAmbiguousAndPartialAreNeverDroppedSilently:
    """The whole point of this module: every input base is accounted for
    in codon_counts, ambiguous_codons, or partial_trailing_bases — never
    just absent from the result."""

    def test_ambiguous_codon_is_counted_not_dropped(self):
        # NNN and one clean codon.
        result = count_codons("ATGNNN")
        assert result.codon_counts == {"ATG": 1}
        assert result.ambiguous_codons == 1
        assert result.counted_codons == 1
        assert result.total_codons == 2  # 1 clean + 1 ambiguous, nothing lost

    def test_partial_trailing_bases_are_reported(self):
        # 7 bases: two full codons + 1 leftover base.
        result = count_codons("ATGATGA")
        assert result.counted_codons == 2
        assert result.partial_trailing_bases == 1

    def test_two_leftover_bases_are_reported(self):
        # 8 bases: two full codons + 2 leftover bases.
        result = count_codons("ATGATGAT")
        assert result.counted_codons == 2
        assert result.partial_trailing_bases == 2

    def test_every_base_is_accounted_for(self):
        seq = "ATGNNNATCAT"  # 11 bases: ATG, NNN, ATC counted/ambiguous, "AT" leftover
        result = count_codons(seq)
        consumed = (result.counted_codons + result.ambiguous_codons) * 3
        assert consumed + result.partial_trailing_bases == len(seq)


class TestRSCU:
    """RSCU compares a codon's use to its synonymous family, and must stay
    silent (not fabricate 0.0) when there's no family data at all."""

    def test_single_codon_amino_acid_is_always_rscu_one(self):
        # ATG (Met) and TGG (Trp) are the only codons for their amino acids,
        # so whenever they're used at all their RSCU is exactly 1.0.
        result = count_codons("ATGTGG")
        assert result.rscu["ATG"] == pytest.approx(1.0)
        assert result.rscu["TGG"] == pytest.approx(1.0)

    def test_even_split_across_a_family_gives_rscu_one(self):
        # Phe (F) has two codons, TTT and TTC. Used equally, both are
        # exactly at the expected rate.
        result = count_codons("TTTTTC")
        assert result.rscu["TTT"] == pytest.approx(1.0)
        assert result.rscu["TTC"] == pytest.approx(1.0)

    def test_skewed_usage_is_reflected_in_rscu(self):
        # Phe family used 4 times, all as TTT: TTT should read above 1.0,
        # its never-used partner TTC should read exactly 0.0 (real data,
        # not missing data — the family total is nonzero).
        result = count_codons("TTTTTTTTTTTT")  # 4x TTT
        assert result.rscu["TTT"] == pytest.approx(2.0)
        assert result.rscu["TTC"] == pytest.approx(0.0)

    def test_unobserved_amino_acid_family_has_no_rscu_entry(self):
        # Only Met (ATG) appears; every other amino acid family has zero
        # observations and must not appear in rscu at all.
        result = count_codons("ATG")
        assert set(result.rscu.keys()) == {"ATG"}
        assert "TTT" not in result.rscu


class TestCodonUsageForOrf:
    def test_orf_codon_usage_starts_at_its_start_codon(self):
        sequence = "ATGAAAAAAAAATGA"
        orfs = find_orfs(sequence, min_length=9, include_reverse=False)
        assert orfs, "fixture ORF must be found for this test to mean anything"

        usage = codon_usage_for_orf(orfs[0])
        assert usage.codon_counts.get("ATG", 0) >= 1
        # The ORF's own sequence ends on a stop codon.
        assert usage.codon_counts.get("TGA", 0) >= 1


class TestMergeCodonUsage:
    def test_merge_sums_across_results(self):
        a = count_codons("ATGATG")
        b = count_codons("ATGATC")
        merged = merge_codon_usage([a, b])

        assert merged.codon_counts["ATG"] == 3
        assert merged.codon_counts["ATC"] == 1
        assert merged.counted_codons == a.counted_codons + b.counted_codons

    def test_merge_recomputes_rscu_over_combined_counts(self):
        # Individually, each half looks 100% one codon; merged, the family
        # is balanced and RSCU should reflect the combined picture, not an
        # average of the two skewed halves.
        a = count_codons("TTTTTT")  # 2x TTT, 0x TTC
        b = count_codons("TTCTTC")  # 0x TTT, 2x TTC
        merged = merge_codon_usage([a, b])

        assert merged.rscu["TTT"] == pytest.approx(1.0)
        assert merged.rscu["TTC"] == pytest.approx(1.0)

    def test_merge_of_empty_list_is_a_valid_empty_result(self):
        merged = merge_codon_usage([])
        assert isinstance(merged, CodonUsageResult)
        assert merged.codon_counts == {}
        assert merged.total_codons == 0


class TestCodonTableRows:
    def test_covers_all_64_codons(self):
        result = count_codons("ATGATC")
        rows = get_codon_table_rows(result)
        assert len(rows) == 64
        assert {r["codon"] for r in rows} == set(
            c for c in rows_all_codons()
        )

    def test_unobserved_family_reports_none_not_zero(self):
        # Only ATG appears, so every other amino acid's rows must show
        # None for fraction/rscu rather than a fabricated 0.
        result = count_codons("ATG")
        rows = get_codon_table_rows(result)

        atg_row = next(r for r in rows if r["codon"] == "ATG")
        assert atg_row["fraction_of_aa_pct"] == pytest.approx(100.0)
        assert atg_row["rscu"] == pytest.approx(1.0)

        untouched_row = next(r for r in rows if r["codon"] == "TTT")
        assert untouched_row["count"] == 0
        assert untouched_row["fraction_of_aa_pct"] is None
        assert untouched_row["rscu"] is None

    def test_rows_sorted_by_amino_acid_then_descending_count(self):
        result = count_codons("TTTTTTTTTTTC")  # 3x TTT, 1x TTC (Phe family)
        rows = get_codon_table_rows(result)
        phe_rows = [r for r in rows if r["amino_acid"] == "F"]
        assert phe_rows[0]["codon"] == "TTT"
        assert phe_rows[0]["count"] >= phe_rows[1]["count"]


def rows_all_codons():
    from app.core.orf import CODON_TABLE
    return CODON_TABLE.keys()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
