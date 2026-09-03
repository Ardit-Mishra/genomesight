"""Coverage for the codon and alignment modules.

Written when the Streamlit implementation was deleted. Its suite had 20 codon and
14 alignment tests, but against a different API (count_codons, codon_usage_for_orf)
that the FastAPI backend does not expose, so they could not simply be moved.
Deleting them without replacement would have shipped less verification than before
-- these cover the same behaviours against the API that actually runs.
"""
import pytest

from app.core.alignment import perform_alignment
from app.core.codons import calculate_codon_usage, calculate_rscu


class TestCodonCounting:
    def test_counts_in_frame_triplets(self):
        assert calculate_codon_usage("ATGGCC") == {"ATG": 1, "GCC": 1}

    def test_lowercase_is_normalised(self):
        assert calculate_codon_usage("atggcc") == {"ATG": 1, "GCC": 1}

    def test_rna_uracil_is_accepted(self):
        assert calculate_codon_usage("AUGGCC") == {"ATG": 1, "GCC": 1}

    def test_trailing_partial_codon_is_ignored(self):
        """5 bases is one codon plus 2 leftover; the leftover is not a codon."""
        assert calculate_codon_usage("ATGGC") == {"ATG": 1}

    def test_stop_codons_are_excluded_from_sense_counts(self):
        """forward_table holds sense codons only; TAA/TAG/TGA are not among them."""
        counts = calculate_codon_usage("ATGTAATAGTGA")
        assert counts == {"ATG": 1}

    def test_ambiguous_codon_is_not_counted_as_a_sense_codon(self):
        assert "ATN" not in calculate_codon_usage("ATGATN")

    def test_empty_sequence(self):
        assert calculate_codon_usage("") == {}


class TestRSCU:
    """RSCU normalises WITHIN a synonymous family; it is not codon frequency.

    A regression here previously shipped to production as count/total, which gives
    0.1538 for every codon of a 13-codon sequence and is not RSCU at all.
    """

    def test_single_codon_amino_acid_is_always_one(self):
        # Methionine has exactly one codon, so its RSCU is 1.0 by definition.
        assert calculate_rscu("ATG")["ATG"] == 1.0

    def test_sole_observed_codon_of_a_four_codon_family(self):
        # Alanine has 4 synonymous codons; observing only GCC gives 1 / (1/4) = 4.0
        assert calculate_rscu("GCC")["GCC"] == 4.0

    def test_evenly_used_family_averages_one(self):
        """Glycine has 4 codons; using two of them equally gives 2.0 each."""
        rscu = calculate_rscu("GGTGGC")
        assert rscu["GGT"] == 2.0
        assert rscu["GGC"] == 2.0

    def test_family_mean_is_one_over_observed_codons(self):
        """The defining property: RSCU averages 1.0 across a family's observed codons."""
        rscu = calculate_rscu("GGTGGCGGAGGG")  # all four glycine codons, once each
        glycine = [rscu[c] for c in ("GGT", "GGC", "GGA", "GGG")]
        assert sum(glycine) / len(glycine) == pytest.approx(1.0)

    def test_is_not_plain_frequency(self):
        """Guards the exact regression that reached production."""
        rscu = calculate_rscu("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        total_codons = sum(calculate_codon_usage("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG").values())
        # count/total would put every value below 1; real RSCU does not.
        assert max(rscu.values()) > 1.0
        assert rscu["GCC"] != pytest.approx(2 / total_codons)


class TestAlignment:
    def test_identical_sequences_are_100_percent(self):
        assert perform_alignment("ACGT", "ACGT", "global")["identity_percentage"] == 100.0

    def test_single_mismatch(self):
        assert perform_alignment("ACGT", "ACTT", "global")["identity_percentage"] == 75.0

    def test_deletion_introduces_a_gap(self):
        result = perform_alignment("ACGT", "AGT", "global")
        assert "-" in result["aligned_seq1"] + result["aligned_seq2"]
        assert result["identity_percentage"] == 75.0

    def test_aligned_sequences_contain_only_sequence_characters(self):
        """The regression guard.

        Biopython's alignment.format() emits coordinate columns --
        'target            0 ACGTACGTAA 10'. Those labels were once consumed as
        sequence, which reported 27.27% identity where the truth was 90.0%.
        """
        result = perform_alignment("ACGTACGTAA", "ACGTCGTAA", "global")
        for key in ("aligned_seq1", "aligned_seq2"):
            assert set(result[key]) <= set("ACGTUN-"), f"{key} contains non-sequence text: {result[key]!r}"
        assert result["identity_percentage"] == 90.0

    def test_aligned_rows_are_the_same_length(self):
        result = perform_alignment("ACGTACGTAA", "ACGTCGTAA", "global")
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"]) == len(result["match_line"])
