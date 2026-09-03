"""
K-mer counting with an optional native (Cython) accelerator.

GenomeSight always runs correctly on a pure-Python k-mer counter. If the
optional Cython extension (app/core/_kmer_accel.pyx, built per
build_kmer_ext.py) has been compiled for the current platform and Python
version, count_kmers() uses it instead and is meaningfully faster on large
inputs -- see benchmarks/benchmark_kmer.py and
benchmarks/kmer_benchmark_results.json for measured numbers.

The two paths are required to return byte-for-byte identical results
(tests/test_kmer_native.py asserts this on both hand-picked and randomized
inputs). Callers must never guess or silently prefer one path -- always use
count_kmers() and read the returned `backend` so the UI can say plainly
which one ran.

This module never masks a real failure with a plausible-looking result:
- If the native extension simply isn't built, that is not an error --
  it's the documented, always-available fallback, and NATIVE_AVAILABLE
  reports it.
- If the native extension IS present but raises while actually counting,
  that's a real bug and is allowed to propagate rather than being caught
  and quietly downgraded to the Python path.
"""
from collections import Counter
from typing import Any, Dict, List, Tuple

try:
    from app.core import _kmer_accel  # type: ignore[attr-defined]
    NATIVE_AVAILABLE = True
    NATIVE_IMPORT_ERROR: str = ""
except ImportError as exc:
    _kmer_accel = None  # type: ignore[assignment]
    NATIVE_AVAILABLE = False
    NATIVE_IMPORT_ERROR = str(exc)

_VALID_BASES = frozenset("ATCG")


def _count_kmers_python(sequences: List[str], k: int) -> Dict[str, int]:
    """
    Reference pure-Python k-mer counter.

    Semantics (kept in lockstep with app/core/_kmer_accel.pyx's
    count_kmers -- see that file's docstring for the exact contract):
    each sequence is scanned independently and upper-cased; a window
    counts only if every base in it is A, T, C, or G.
    """
    counts: Counter = Counter()
    for sequence in sequences:
        seq = sequence.upper()
        n = len(seq)
        for i in range(n - k + 1):
            kmer = seq[i:i + k]
            if all(base in _VALID_BASES for base in kmer):
                counts[kmer] += 1
    return dict(counts)


def count_kmers(sequences: List[str], k: int) -> Tuple[Dict[str, int], str]:
    """
    Count k-mers of length `k` across `sequences`.

    Args:
        sequences: DNA/RNA sequence strings (already extracted from any
            Bio.SeqRecord wrapper -- this function is Streamlit- and
            Biopython-free).
        k: k-mer length, must be >= 1.

    Returns:
        (counts, backend) where counts maps each observed k-mer to its
        count, and backend is "native" or "python" -- the caller (e.g.
        SequenceAnalyzer.analyze_kmers) is expected to surface `backend`
        to the UI rather than swallow it.

    Raises:
        ValueError: if k < 1.
        RuntimeError: if the native extension is present but fails while
            counting (a real bug -- not silently downgraded to Python).
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")

    if NATIVE_AVAILABLE:
        try:
            counts = _kmer_accel.count_kmers(list(sequences), k)
        except ValueError:
            raise
        except Exception as exc:  # pragma: no cover - defends a real extension bug
            raise RuntimeError(
                f"native k-mer accelerator raised while counting: {exc}"
            ) from exc
        return dict(counts), "native"

    return _count_kmers_python(list(sequences), k), "python"


def backend_status() -> Dict[str, Any]:
    """Machine-readable status for diagnostics / the UI announcement."""
    return {
        "native_available": NATIVE_AVAILABLE,
        "backend": "native" if NATIVE_AVAILABLE else "python",
        "import_error": NATIVE_IMPORT_ERROR,
    }
