# What was decided, what it cost, what is still unknown

The most consequential decision here was to put **no model anywhere in the
factual path**. Everything else follows from it.

---

## Why there is no machine learning in this tool

Every other project in this portfolio trains something. This one deliberately
does not. The operations it performs — GC and base composition, k-mer profiling,
open-reading-frame detection, IUPAC motif and restriction-site search, codon
usage, pairwise alignment — all have exact answers. A model could only introduce
variance into questions that have none.

So the guarantee this tool offers is one the others cannot: **the same input
always returns the same result**, and every result names the method that
produced it rather than presenting a bare number.

An optional classifier over public benchmark data was scoped and explicitly
deferred, on the grounds that it would have to be labelled as a statistical
classifier on benchmark data and could never be the primary move here.

## The Cython accelerator, measured rather than asserted

k-mer counting is the hot path. A native extension was written for it, and the
claim is measured across sequence lengths and k values rather than quoted from a
single flattering run:

| Sequence length | k | Python (s) | Native (s) | Speed-up |
|---:|---:|---:|---:|---:|
| 1,000 | 4 | 0.000784 | 0.000121 | 6.5× |
| 1,000 | 12 | 0.001230 | 0.000111 | **11.0×** |
| 10,000 | 8 | 0.011464 | 0.001668 | 6.9× |
| 10,000 | 12 | 0.013459 | 0.001226 | 11.0× |
| 100,000 | 8 | 0.138141 | 0.020162 | 6.9× |
| 500,000 | 8 | 0.707789 | 0.109870 | 6.4× |
| 500,000 | 12 | 0.921982 | 0.174901 | **5.3×** |

Source: `backend/benchmarks/kmer_benchmark_results.json`, produced by
`benchmark_kmer.py`. Every row was checked to return the identical unique-k-mer
count from both paths — a faster wrong answer is not an optimisation.

Two things worth reading off that table. The headline is the range, **5.3×–11.0×**,
not the best number in it. And the advantage *shrinks* as input grows, because
both paths still build the same dictionary and that allocation comes to dominate.
Quoting 11× alone would describe the small-input case and mislead about the large
one.

**The pure-Python path is never removed and never silently substituted.** The
engine actually serving a request is reported at runtime, so a failed native
build costs speed, not correctness, and the reader can tell which one ran.

## A bug that changed what "verified" means here

The pairwise alignment module was parsing Biopython's **coordinate columns as
sequence**. It reported **81.1% identity where the true value was 92.9%** — a
plausible-looking number, wrong by twelve points, on a function whose whole
purpose is to be exact.

Nothing about the output looked broken. That is the point: a deterministic tool
returning a confidently wrong number is worse than a model that admits
uncertainty, because there is no uncertainty on display to warn anyone.

The fix came with a regression guard pinning the known-correct identity, so the
same class of parse error fails the suite rather than the user.

There was a related structural problem: a correct alignment implementation
existed in the backend but was **unreachable from the interface**. Working code
nobody can invoke is not a feature, and wiring it through was treated as the
higher-value move than adding anything new.

## What is still unknown

- The benchmark was run on one machine (Windows 11, Python 3.12). The speed-up
  range should be expected to move on other platforms and compilers.
- Alignment is pairwise only; no multiple-sequence alignment.
- Large-input behaviour is measured to 500,000 bases. Beyond that, the dictionary
  allocation that already dominates at that size will dominate further, and no
  measurement backs any claim past it.

## Scope

Research and educational tooling for routine sequence work. It computes what it
says it computes and makes no claim beyond it.
