![Project Banner](assets/banner.png)
# GenomeSight

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A sequence analysis workbench: a FastAPI service and a React front end for exploratory
analysis of DNA and RNA. Composition and GC content, k-mer profiling, six-frame ORF
detection, IUPAC motif and restriction-site search, codon usage with RSCU, translation,
and pairwise alignment.

**Live:** https://genomesight-frontend.vercel.app · **API:** https://genomesight-api.onrender.com

No model sits in the factual path. The same input always produces the same result.

## Why this exists

Routine sequence work — checking GC content, finding open reading frames, looking for a
restriction site, comparing two reads — usually means either writing a throwaway Biopython
script or pasting sequence into a site that offers no account of how it computed anything.

GenomeSight does these operations in one place and makes each result checkable: the method
that ran is named, the alignment identity is defined over aligned columns, the k-mer engine
reports whether the native or Python path produced the counts, and the API states the commit
it was built from. It is aimed at anyone who wants a quick answer they can verify rather than
one they must trust.

## Features

- **Sequence statistics** — GC content, nucleotide composition, molecular weight; FASTA,
  FASTQ and GenBank input, detected by content rather than by file extension
- **K-mer analysis** — configurable k, with an optional Cython accelerator measured at
  **5.27×–11.04×** over the pure-Python path (see [Native k-mer accelerator](#native-k-mer-accelerator))
- **ORF detection** — all three forward frames plus the reverse complement, with a
  minimum-length threshold and per-frame summary
- **Motif and restriction-site search** — IUPAC ambiguity codes expanded to a regular
  expression that the API returns, so the expansion is checkable rather than trusted
- **Codon usage** — counts and relative synonymous codon usage, normalised within each
  amino acid's synonymous family
- **Pairwise alignment** — global Needleman–Wunsch with identity over aligned columns
- **Export** — every results panel exports to CSV, containing all rows rather than the
  subset displayed

## Repository layout

```
backend/          FastAPI service
  app/core/       analysis modules (no web framework imports)
  app/api/        routes and schemas
  benchmarks/     k-mer benchmark and its recorded results
  tests/          67 tests
  Dockerfile      builds the native accelerator, honours $PORT
frontend/         React + Vite + Tailwind
sample_data/      example FASTA/FASTQ files
DEPLOYMENT.md     hosts, health check, rollback, incident notes
```

## Running locally

**Backend** (Python 3.12):

```bash
cd backend
uv run --python 3.12 --with-requirements requirements.txt uvicorn app.main:app --reload --port 8000
```

**Frontend** (Node 20+):

```bash
cd frontend
npm install
echo "VITE_API_URL=http://127.0.0.1:8000" > .env.local
npm run dev
```

Interactive API docs are at `http://localhost:8000/docs`.

## API

| Method | Path | Purpose |
|---|---|---|
| POST | `/api/analyze` | composition, GC, k-mer profile |
| POST | `/api/align` | pairwise global alignment |
| POST | `/api/translate` | frame translation |
| POST | `/api/codons` | codon counts and RSCU |
| POST | `/api/orfs` | six-frame ORF detection |
| POST | `/api/motifs` | IUPAC pattern search |
| POST | `/api/restriction-sites` | named restriction enzyme |
| GET | `/api/enzymes` | known recognition sites |
| GET | `/api/health` | liveness, **deployed commit**, k-mer engine |

`/api/health` reports the commit the running image was built from. That field exists
because a deployment once ran five commits behind for hours and nothing served by the
API could distinguish a stale deploy from a code defect. See DEPLOYMENT.md.

## Native k-mer accelerator

K-mer counting has two paths:

1. **Pure Python** — `backend/app/core/kmer_native._count_kmers_python`, always available
2. **Native Cython** — `backend/app/core/_kmer_accel.pyx`, optional

`count_kmers()` uses the compiled extension when present and the Python path otherwise,
and returns which one ran. The API surfaces that as `kmer_engine`, so the active path is
always visible rather than assumed. The two are required to return identical results, and
a test asserts it on hand-picked and randomised inputs including ambiguous bases.

Build it (needs a C compiler, Cython and setuptools — none required at runtime):

```bash
cd backend
uv run --python 3.12 --with cython --with setuptools --with-requirements requirements.txt \
    python build_kmer_ext.py build_ext --inplace
```

The Docker image builds it automatically. A failed build does not fail the image: the
service runs the Python path and says so.

### Benchmark

```bash
cd backend
uv run --python 3.12 --with-requirements requirements.txt python benchmarks/benchmark_kmer.py
```

Measured on one Windows desktop across 12 configurations: **5.27×–11.04×**. The speedup
shrinks as input grows, because both paths still build the same Python dictionary — the
accelerator removes per-window interpreter overhead, not the dictionary. Recorded in
`backend/benchmarks/kmer_benchmark_results.json`.

## Limitations

Stated plainly, because a tool that claims no weaknesses invites the reader to find them.

- **Alignment is pairwise only.** No multiple sequence alignment. The implementation wraps
  Biopython's `PairwiseAligner` with default scoring, which is not tuned for any particular
  biological question.
- **ORF coordinates on the minus strand** are reported against the reverse-complement
  sequence, not the original strand. The API returns a `coordinate_note` saying so, but a
  caller that ignores it will misplace features.
- **ORF detection is codon-based, not gene prediction.** It finds start-to-stop spans. It
  does not model splicing, ribosomal binding sites, or coding potential, so it will report
  spans that are not genes.
- **`find_common_motifs` counts exact k-mers.** It does not build a position weight matrix
  and does not score statistical enrichment against a background model.
- **The restriction enzyme set is 10 common enzymes**, hardcoded. It is not REBASE.
- **The k-mer benchmark is single-machine, single-run.** 12 configurations on one Windows
  desktop, best-of-N timing, no confidence intervals. It shows the native path is faster
  here; it is not a portable performance claim.
- **The hosted backend sleeps.** On the free tier the first request after ~15 minutes idle
  waits ~50 seconds. The client retries transparently, but the wait is real.
- **Large inputs are held in memory.** There is no streaming parser, so genome-scale FASTA
  will exhaust memory rather than degrade gracefully.

## Reproduce

Every number in this README comes from a committed artifact:

| claim | regenerate with |
|---|---|
| 67 tests | `cd backend && pytest tests/ -q` |
| 5.27×–11.04× k-mer speedup | `cd backend && python benchmarks/benchmark_kmer.py` |
| the deployed commit | `curl https://genomesight-api.onrender.com/api/health` |

## Tests

```bash
cd backend
uv run --python 3.12 --with-requirements requirements.txt --with-requirements requirements-dev.txt \
    -m pytest tests/ -q
```

67 tests. `httpx` is pinned below 0.28 in `requirements-dev.txt`: Starlette's `TestClient`
calls `httpx.Client(app=...)`, removed in 0.28, and the result is a **collection** error —
the suite does not fail, it fails to run, reporting zero tests.

Several tests are regression guards for defects that reached production:

- Alignment once parsed Biopython's coordinate columns as sequence, reporting 81.1%
  identity where the truth was 92.9%, and 27.27% where it was 90.0%
- RSCU once returned `count / total_codons`, which is codon frequency, not RSCU
- The k-mer table was captioned "ranked by count" while rendering insertion order

## Acknowledgements

- [Biopython](https://biopython.org/) for sequence parsing and alignment
- [FastAPI](https://fastapi.tiangolo.com/) and [React](https://react.dev/)
- [Cython](https://cython.org/) for the optional accelerator

## License

MIT — see [LICENSE](LICENSE).
