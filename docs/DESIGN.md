# GenomeSight — architecture and design decisions

**Original design:** 2026-09-02 · **Reconciled against the build:** 2026-09-03

This records why GenomeSight is shaped the way it is. Where the built system diverged from
the original design, both are shown — a design document that quietly matches whatever was
built teaches nothing, and one that silently doesn't match misleads.

---

## 1. Objectives

GenomeSight moved from a single Streamlit script to a decoupled full-stack application: a
stateless FastAPI service holding all the Biopython computation, and a React single-page
frontend that talks to it over REST.

Three shifts drove it:

1. **Frontend/backend decoupling.** The analysis code has no web framework imports, so it
   is testable directly and callable by anything — not only by a UI.
2. **Cold-start mitigation.** The backend runs on a free tier that sleeps when idle. The
   client retries with exponential backoff rather than surfacing a failure, and an external
   scheduled job keeps the service warm.
3. **A scientific-instrument interface** rather than a dashboard: monospace sequence data,
   restrained colour, and results that state their own provenance.

## 2. Structure

```
backend/
  app/
    main.py              FastAPI app, CORS, service index
    api/
      endpoints.py       routes
      schemas.py         Pydantic request/response models
    core/                analysis modules — no web framework imports
      alignment.py       pairwise global/local alignment
      io_fasta.py        FASTA / FASTQ / GenBank parsing, format detected by content
      codons.py          codon usage and RSCU
      kmer_native.py     k-mer counting, native and pure-Python paths
      _kmer_accel.pyx    optional Cython accelerator
      orf.py             six-frame ORF detection
      motifs.py          IUPAC motif and restriction-site search
      validation.py      nucleotide validation and normalisation
      build_info.py      the commit the running image was built from
      pinger.py          process liveness
  benchmarks/            k-mer benchmark and recorded results
  tests/                 67 tests
frontend/                React + Vite + Tailwind
```

**Diverged from the original design.** The design named a `sequence_analyzer.py` holding
composition, entropy and stats. It was never built: those computations turned out small
enough to sit in the endpoint alongside the k-mer call, and a module existing only to group
three functions is indirection without benefit. `orf.py`, `motifs.py` and `build_info.py`
post-date the design.

## 3. API

| endpoint | purpose |
|---|---|
| `GET /` | service index, enumerated from the router so it cannot drift |
| `GET /api/health` | liveness, **deployed commit**, active k-mer engine |
| `POST /api/analyze` | composition, GC, k-mer profile |
| `POST /api/align` | pairwise alignment, identity over aligned columns |
| `POST /api/translate` | frame translation |
| `POST /api/codons` | codon usage and RSCU |
| `POST /api/orfs` | six-frame ORF detection |
| `POST /api/motifs` | IUPAC pattern search |
| `POST /api/restriction-sites` | named restriction enzyme |
| `GET /api/enzymes` | known recognition sites |

**`/api/health` reporting the deployed commit was not designed.** It was added after a
deployment ran five commits behind for hours while every symptom pointed at the code being
wrong. Nothing the API served identified which commit it was built from, so telling a stale
deploy from a defect meant comparing behaviour endpoint by endpoint. One field replaced
that entirely.

## 4. Decisions worth stating

**No model sits in the factual path.** Every result is computed deterministically, so the
same input always produces the same output and any result can be checked by hand.

**The k-mer accelerator reports which path ran.** There are two implementations — a Cython
extension and a pure-Python fallback — required to return identical results, asserted by a
test on randomised input including ambiguous bases. The API returns `kmer_engine`, so the
active path is visible rather than assumed. A silent fallback turns a measured performance
claim into an unfalsifiable one.

**Format detection reads content, not filenames.** A GenBank file begins with `LOCUS`;
testing only for `>` and `@` silently rejects it, which is what an earlier version did while
the documentation claimed GenBank support.

**Alignment identity is computed over aligned columns only.** Biopython's `format()` emits
coordinate columns — `target  0 ACGTACGTAA 10`. Consuming those as sequence reported 27.27%
identity where the truth was 90.0%. A regression test now asserts the aligned strings
contain nothing but sequence characters.

## 5. Deployment

- **Frontend:** Vercel static build.
- **Backend:** Render free container, deploying from `main`.
- **Cold start:** the instance sleeps after ~15 minutes idle and the next request waits
  ~50 seconds. The client retries with exponential backoff, and a scheduled external job
  pings `/api/health` every 12 minutes. A self-ping cannot solve this — once the instance
  sleeps there is no process left to send one.

See [DEPLOYMENT.md](../DEPLOYMENT.md) for hosts, rollback, and the branch-mismatch incident.
