# cython: language_level=3, boundscheck=False, wraparound=False
"""
Cython-accelerated k-mer counter for GenomeSight.

This module is a compiled C extension, NOT part of the app's required
dependency set. It is built explicitly (see build_kmer_ext.py at the repo
root and README.md's "Native k-mer accelerator" section) and is entirely
optional: app/core/kmer_native.py imports it if present and falls back to
a pure-Python implementation with byte-for-byte identical semantics when
it is not built for the current platform/Python version.

Semantics (must match app.core.kmer_native._count_kmers_python exactly):
    - Each sequence in `sequences` is upper-cased and scanned independently
      (k-mers never span a boundary between two sequences).
    - A window of length k is counted only if every base in it is one of
      A, T, C, G (anything else -- N, IUPAC ambiguity codes, whitespace,
      lowercase-after-upper artifacts, non-ASCII -- invalidates that window,
      exactly as the pure-Python path does).
    - Counting is unordered (a plain dict of kmer -> count).
"""

from libc.stdint cimport uint8_t


cpdef dict count_kmers(list sequences, int k):
    """
    Count valid ATCG k-mers of length `k` across `sequences`.

    Args:
        sequences: list of DNA/RNA sequence strings.
        k: k-mer length, must be >= 1.

    Returns:
        dict mapping each observed k-mer (str) to its occurrence count (int).

    Raises:
        ValueError: if k < 1.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")

    cdef dict counts = {}
    cdef object sequence
    cdef bytes seq_bytes
    cdef const unsigned char[:] view
    cdef Py_ssize_t n, i, j
    cdef bint valid
    cdef unsigned char base
    cdef str kmer

    for sequence in sequences:
        # encode with 'replace' (not 'ignore') so every character -- ASCII
        # or not -- still occupies one byte. That keeps window positions
        # identical to the pure-Python path, which indexes the original
        # (unencoded) string one character at a time.
        seq_bytes = (<str> sequence).upper().encode("ascii", "replace")
        n = len(seq_bytes)
        if n < k:
            continue
        view = seq_bytes

        for i in range(n - k + 1):
            valid = True
            for j in range(i, i + k):
                base = view[j]
                if base != 65 and base != 84 and base != 67 and base != 71:
                    # 65='A' 84='T' 67='C' 71='G'
                    valid = False
                    break
            if valid:
                kmer = seq_bytes[i:i + k].decode("ascii")
                if kmer in counts:
                    counts[kmer] += 1
                else:
                    counts[kmer] = 1

    return counts
