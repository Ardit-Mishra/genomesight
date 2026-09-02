# GenomeSight FastAPI + React Rewrite — Design Specification

**Date:** 2026-09-02  
**Status:** Approved  
**Project:** GenomeSight (`genomesight.arditmishra.com`)  
**Stack:** Python FastAPI (Backend) + React / Vite / Tailwind CSS / Lucide Icons / Chart.js (Frontend) + Render & Vercel ($0 Deployment)

---

## 1. Executive Summary & Objectives

GenomeSight is being upgraded from a prototype Streamlit script to a production-grade, decoupled full-stack architecture. This showcases professional software engineering, full-stack separation of concerns, and clean API design for recruiter evaluation (AI/ML & Bioinformatics roles).

### Key Architectural Shifts:
1. **Frontend/Backend Decoupling:** Stateless REST API via FastAPI handling all Biopython computation, decoupled entirely from the UI. Single-page React application hosted on Vercel.
2. **Eliminating Cold-Start Friction:** Render backend paired with a lightweight background pinger / cron heartbeat + resilient client-side retry handling.
3. **Industrial Utilitarian UI/UX:** Dark-mode genomic laboratory theme with high-contrast typography (IBM Plex Mono + Inter), monospace data readouts, and responsive interactive visualizations.

---

## 2. Monorepo & Folder Structure

```
genomesight/
├── backend/
│   ├── app/
│   │   ├── __init__.py
│   │   ├── main.py              # FastAPI app initialization, CORS, health check
│   │   ├── api/
│   │   │   ├── __init__.py
│   │   │   ├── endpoints.py     # API routers (/analyze, /alignment, /codons, /translate)
│   │   │   └── schemas.py       # Pydantic request/response models
│   │   └── core/
│   │       ├── __init__.py
│   │       ├── sequence_analyzer.py # GC content, composition, entropy, stats
│   │       ├── alignment.py     # Needleman-Wunsch & Smith-Waterman wrappers
│   │       ├── io_fasta.py      # FASTA/FASTQ/GenBank parsers via BioPython
│   │       ├── codons.py        # Codon usage and RSCU calculator
│   │       ├── kmer_native.py   # K-mer frequency analysis
│   │       └── pinger.py        # Self-ping / warm-up utilities
│   ├── requirements.txt
│   └── Dockerfile               # For Render container deployment
├── frontend/
│   ├── public/
│   ├── src/
│   │   ├── components/
│   │   │   ├── Navbar.tsx
│   │   │   ├── FileUpload.tsx
│   │   │   ├── SummaryCards.tsx
│   │   │   ├── CompositionCharts.tsx
│   │   │   ├── AlignmentViewer.tsx
│   │   │   └── CodonTable.tsx
│   │   ├── services/
│   │   │   └── api.ts           # Axios client with retry logic for Render cold start
│   │   ├── App.tsx
│   │   ├── index.css
│   │   └── main.tsx
│   ├── package.json
│   ├── tailwind.config.js
│   └── vite.config.ts
├── docs/
│   └── superpowers/
│       └── specs/
│           └── 2026-09-02-genomesight-fastapi-react-design.md
├── DEPLOYMENT.md
└── README.md
```

---

## 3. Backend API Specifications (FastAPI)

### Endpoints:
- `GET /api/health` — Health check returning status and backend version.
- `POST /api/analyze` — Accepts FASTA text or file upload; computes sequence composition, lengths, GC percentages, complexity entropy, and k-mer distributions.
- `POST /api/align` — Accepts `seq1`, `seq2`, and mode (`global` or `local`), returning gapped alignments, match lines, identity percentage, and score.
- `POST /api/translate` — Translates nucleotide sequences to amino acid chains.
- `POST /api/codons` — Computes Codon usage and Relative Synonymous Codon Usage (RSCU).

---

## 4. Frontend Design & UI/UX Standards

- **Theme:** Industrial Utilitarian / Genomic Laboratory (Dark mode `#0a0a0f` canvas, `#14141d` surface, `#22c55e` bio-green accent, `#38bdf8` cyan telemetry accents).
- **Typography:** `IBM Plex Mono` for sequence data and metrics; `Inter` for UI labels.
- **Components:**
  - Drag-and-drop FASTA file uploader with sample sequence quick-load.
  - Summary cards with key metrics (total length, average GC%, N50/entropy).
  - Dual-mode visualization: Rich data tables + interactive charts (Chart.js) for composition + Zoomable SVG alignment viewer for sequence matching.

---

## 5. Deployment Strategy ($0 Marginal Cost)

- **Frontend:** Vercel static deployment (`genomesight.arditmishra.com`).
- **Backend:** Render free container web service (`genomesight-api.onrender.com`).
- **Cold-Start Mitigation:** Client-side Axios interceptor with automatic retry (up to 3 attempts with exponential backoff and a "Warming up bioinformatics engine..." toast notification).
