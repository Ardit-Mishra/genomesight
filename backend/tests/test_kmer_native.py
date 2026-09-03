"""
Tests for the k-mer counting dispatcher (app.core.kmer_native) and, when
built, the native Cython accelerator (app.core._kmer_accel).

Run with: pytest tests/test_kmer_native.py -v

These tests are the honesty backstop for the native-accelerator claim:
if the compiled extension and the pure-Python counter ever disagree, that
is a real bug and these tests must fail loudly, never pass silently.
"""
import random

import pytest

from app.core import kmer_native


class TestPythonCounterBaseline:
    """The pure-Python path must always work -- it's the required fallback."""

    def test_counts_simple_sequence(self):
        counts = kmer_native._count_kmers_python(["ATGCATGC"], 3)
        assert counts == {
            "ATG": 2, "TGC": 2, "GCA": 1, "CAT": 1,
        }

    def test_case_insensitive(self):
        lower = kmer_native._count_kmers_python(["atgcatgc"], 3)
        upper = kmer_native._count_kmers_python(["ATGCATGC"], 3)
        assert lower == upper

    def test_ambiguous_base_breaks_windows(self):
        # "N" is not A/T/C/G, so every 3-mer touching it is invalid.
        counts = kmer_native._count_kmers_python(["ATNGCA"], 3)
        assert counts == {"GCA": 1}

    def test_sequences_do_not_bleed_into_each_other(self):
        # A k-mer must never span the boundary between two input sequences.
        joined = kmer_native._count_kmers_python(["ATG"], 3)
        separate = kmer_native._count_kmers_python(["AT", "G"], 3)
        assert joined == {"ATG": 1}
        assert separate == {}

    def test_short_sequence_yields_no_kmers(self):
        assert kmer_native._count_kmers_python(["AT"], 3) == {}

    def test_empty_sequence_list(self):
        assert kmer_native._count_kmers_python([], 3) == {}


class TestCountKmersDispatcher:
    """count_kmers() must report which backend actually ran."""

    def test_backend_matches_availability(self):
        counts, backend = kmer_native.count_kmers(["ATGCATGC"], 3)
        if kmer_native.NATIVE_AVAILABLE:
            assert backend == "native"
        else:
            assert backend == "python"
        assert counts == {"ATG": 2, "TGC": 2, "GCA": 1, "CAT": 1}

    def test_invalid_k_raises(self):
        with pytest.raises(ValueError):
            kmer_native.count_kmers(["ATGC"], 0)

    def test_backend_status_reports_reality(self):
        status = kmer_native.backend_status()
        assert status["native_available"] == kmer_native.NATIVE_AVAILABLE
        assert status["backend"] in ("native", "python")
        if not kmer_native.NATIVE_AVAILABLE:
            assert status["import_error"], (
                "native accelerator absent but no import error was recorded -- "
                "the UI has nothing honest to report"
            )


@pytest.mark.skipif(
    not kmer_native.NATIVE_AVAILABLE,
    reason=(
        "native k-mer accelerator not built in this environment "
        f"({kmer_native.NATIVE_IMPORT_ERROR}) -- run "
        "`python build_kmer_ext.py build_ext --inplace` (see README) to enable it"
    ),
)
class TestNativeMatchesPython:
    """
    The whole point of shipping two implementations is that they must
    agree. Any mismatch here means the native path is wrong -- not that
    it should be silently preferred or silently ignored.
    """

    def _random_sequence(self, rng: random.Random, length: int, alphabet: str) -> str:
        return "".join(rng.choice(alphabet) for _ in range(length))

    def test_agrees_on_hand_picked_cases(self):
        cases = [
            (["ATGCATGC"], 3),
            (["atgcATGCatgc"], 4),
            (["ATNGCAN"], 2),
            (["A"], 1),
            (["AAAAAAAAAA"], 5),
            (["ATG", "CGT", "AAAT"], 2),
            (["ACGTACGTACGTN"], 6),
        ]
        for sequences, k in cases:
            python_counts = kmer_native._count_kmers_python(sequences, k)
            native_counts = kmer_native._kmer_accel.count_kmers(list(sequences), k)
            assert native_counts == python_counts, (sequences, k)

    def test_agrees_on_randomized_inputs(self):
        rng = random.Random(20260830)
        alphabets = ["ATCG", "ATCGN", "ATCGNRYSWKM", "atcgATCG"]
        for _ in range(30):
            num_seqs = rng.randint(1, 4)
            k = rng.randint(1, 8)
            sequences = [
                self._random_sequence(rng, rng.randint(0, 40), rng.choice(alphabets))
                for _ in range(num_seqs)
            ]
            python_counts = kmer_native._count_kmers_python(sequences, k)
            native_counts = kmer_native._kmer_accel.count_kmers(list(sequences), k)
            assert native_counts == python_counts, (sequences, k)

    def test_native_invalid_k_raises(self):
        with pytest.raises(ValueError):
            kmer_native._kmer_accel.count_kmers(["ATGC"], 0)

    def test_dispatcher_actually_uses_native(self):
        counts, backend = kmer_native.count_kmers(["ATGCATGC"], 3)
        assert backend == "native"
        assert counts == kmer_native._count_kmers_python(["ATGCATGC"], 3)
