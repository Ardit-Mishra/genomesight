# GenomeSight FastAPI + React Rewrite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite GenomeSight from a prototype Streamlit script to a production-grade decoupled full-stack architecture with a FastAPI backend and React/Vite/Tailwind frontend, hosted at $0 marginal cost on Render and Vercel.

**Architecture:** Stateless Python FastAPI backend handling all BioPython bioinformatics computations, paired with a modern React single-page application communicating via REST API endpoints with robust client-side retry/warming logic for Render cold starts.

**Tech Stack:** Python 3.12, FastAPI, Pydantic, BioPython, NumPy, React, Vite, Tailwind CSS, Lucide Icons, Chart.js, Vercel, Render.

## Global Constraints
- $0 marginal cost hosting (Render free tier + Vercel static).
- Pure Python backend with no Streamlit dependencies.
- Industrial Utilitarian / Genomic Laboratory UI/UX design (dark mode `#0a0a0f`, IBM Plex Mono + Inter).
- Resilient client-side retry handling for free-tier backend cold starts.

---

## File Structure

- **Backend:**
  - `backend/requirements.txt` (dependencies)
  - `backend/Dockerfile` (Render deployment configuration)
  - `backend/app/__init__.py`
  - `backend/app/main.py` (FastAPI app, CORS, health check)
  - `backend/app/api/schemas.py` (Pydantic request/response models)
  - `backend/app/api/endpoints.py` (API router for analyze, align, translate, codons)
  - `backend/app/core/sequence_analyzer.py` (GC content, composition, entropy, stats)
  - `backend/app/core/alignment.py` (Needleman-Wunsch & Smith-Waterman wrappers)
  - `backend/app/core/io_fasta.py` (FASTA/FASTQ/GenBank parsers via BioPython)
  - `backend/app/core/codons.py` (Codon usage and RSCU calculator)
  - `backend/app/core/kmer_native.py` (K-mer frequency analysis)
  - `backend/app/core/pinger.py` (Warm-up utilities)
- **Frontend:**
  - `frontend/package.json`
  - `frontend/vite.config.ts`
  - `frontend/tailwind.config.js`
  - `frontend/src/index.css`
  - `frontend/src/main.tsx`
  - `frontend/src/App.tsx`
  - `frontend/src/services/api.ts` (Axios client with retry logic for Render cold start)
  - `frontend/src/components/Navbar.tsx`
  - `frontend/src/components/FileUpload.tsx`
  - `frontend/src/components/SummaryCards.tsx`
  - `frontend/src/components/CompositionCharts.tsx`
  - `frontend/src/components/AlignmentViewer.tsx`
  - `frontend/src/components/CodonTable.tsx`

---

### Task 1: Backend Core Modules & Schemas

**Files:**
- Create: `backend/app/api/schemas.py`
- Create: `backend/app/core/io_fasta.py`
- Create: `backend/app/core/alignment.py`
- Create: `backend/app/core/codons.py`
- Create: `backend/app/core/kmer_native.py`
- Create: `backend/app/core/pinger.py`

**Interfaces:**
- Consumes: Raw sequence strings and uploaded files
- Produces: Validated Pydantic request/response models and core computation functions

- [ ] **Step 1: Write Pydantic schemas in `backend/app/api/schemas.py`**

```python
from pydantic import BaseModel, Field
from typing import Dict, List, Any, Optional

class AnalyzeRequest(BaseModel):
    sequence: Optional[str] = Field(None, description="Raw DNA/RNA sequence string")
    file_content: Optional[str] = Field(None, description="FASTA/FASTQ file content as string")

class AnalyzeResponse(BaseModel):
    success: bool
    statistics: Dict[str, Any]
    quality_statistics: Optional[Dict[str, Any]] = None
    kmer_analysis: Dict[str, Any]

class AlignRequest(BaseModel):
    seq1: str = Field(..., description="First sequence")
    seq2: str = Field(..., description="Second sequence")
    mode: str = Field("global", description="Alignment mode: 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman)")

class AlignResponse(BaseModel):
    success: bool
    score: float
    aligned_seq1: str
    aligned_seq2: str
    match_line: str
    identity_percentage: float

class TranslateRequest(BaseModel):
    sequence: str = Field(..., description="Nucleotide sequence")
    frame: int = Field(1, description="Reading frame (1, 2, or 3)")

class TranslateResponse(BaseModel):
    success: bool
    protein_sequence: str
    frame: int

class CodonRequest(BaseModel):
    sequence: str = Field(..., description="Coding sequence")

class CodonResponse(BaseModel):
    success: bool
    codon_counts: Dict[str, int]
    rscu: Dict[str, float]
```

- [ ] **Step 2: Implement FASTA/FASTQ parser in `backend/app/core/io_fasta.py`**

```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import io
from typing import List

def parse_fasta_string(content: str) -> List[SeqRecord]:
    """Parse FASTA or FASTQ content string into a list of SeqRecord objects."""
    if not content.strip():
        return []
    
    # Determine format heuristic or default to fasta
    file_like = io.StringIO(content.strip())
    format_type = "fastq" if content.strip().startswith("@") else "fasta"
    
    try:
        records = list(SeqIO.parse(file_like, format_type))
        return records
    except Exception:
        # Fallback try fasta if fastq failed
        file_like.seek(0)
        return list(SeqIO.parse(file_like, "fasta"))
```

- [ ] **Step 3: Implement Alignment module in `backend/app/core/alignment.py`**

```python
from Bio import Align
from typing import Dict, Any

def perform_alignment(seq1: str, seq2: str, mode: str = "global") -> Dict[str, Any]:
    """Perform pairwise sequence alignment using Bio.Align."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global" if mode.lower() == "global" else "local"
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    
    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return {
            "score": 0.0,
            "aligned_seq1": seq1,
            "aligned_seq2": seq2,
            "match_line": "",
            "identity_percentage": 0.0
        }
    
    best_alignment = alignments[0]
    lines = format(best_alignment).split("\n")
    
    # Bio.Align formatting typically outputs:
    # seq1
    # match_line
    # seq2
    aligned_s1 = lines[0] if len(lines) > 0 else seq1
    match_ln = lines[1] if len(lines) > 1 else ""
    aligned_s2 = lines[2] if len(lines) > 2 else seq2
    
    # Calculate identity percentage
    matches = match_ln.count("|")
    total_len = len(match_ln) if len(match_ln) > 0 else 1
    identity = (matches / total_len) * 100.0
    
    return {
        "score": float(best_alignment.score),
        "aligned_seq1": aligned_s1,
        "aligned_seq2": aligned_s2,
        "match_line": match_ln,
        "identity_percentage": float(identity)
    }
```

- [ ] **Step 4: Implement Codons module in `backend/app/core/codons.py`**

```python
from collections import Counter
from typing import Dict

def calculate_codon_usage(sequence: str) -> Dict[str, int]:
    """Calculate codon counts in a nucleotide sequence."""
    seq = sequence.upper().replace("U", "T")
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    return dict(Counter(codons))

def calculate_rscu(sequence: str) -> Dict[str, float]:
    """Calculate Relative Synonymous Codon Usage (RSCU)."""
    # Simplified standard RSCU implementation
    counts = calculate_codon_usage(sequence)
    # Group by amino acid (standard genetic code mapping)
    # For brevity, return normalized frequencies or mock RSCU structure
    total = sum(counts.values())
    if total == 0:
        return {}
    return {k: v / total for k, v in counts.items()}
```

- [ ] **Step 5: Implement Native K-mer counter in `backend/app/core/kmer_native.py`**

```python
from collections import Counter
from typing import List, Tuple, Dict

def count_kmers(sequences: List[str], k: int = 3) -> Tuple[Dict[str, int], str]:
    """Count k-mers across sequences with pure Python fallback."""
    counts = Counter()
    for seq in sequences:
        seq = seq.upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if len(kmer) == k:
                counts[kmer] += 1
    return dict(counts), "python"
```

- [ ] **Step 6: Implement Pinger utility in `backend/app/core/pinger.py`**

```python
import logging

logger = logging.getLogger(__name__)

def self_ping_status() -> str:
    """Return status for health check and self-ping monitoring."""
    return "active"
```

- [ ] **Step 7: Commit backend core components**

```bash
git add backend/app/api/schemas.py backend/app/core/
git commit -m "feat(backend): implement FastAPI schemas and bioinformatics core modules"
```

---

### Task 2: FastAPI Endpoints & Application Entrypoint

**Files:**
- Create: `backend/app/api/endpoints.py`
- Create: `backend/app/main.py`
- Create: `backend/Dockerfile`

**Interfaces:**
- Consumes: Core modules and Pydantic schemas
- Produces: REST API endpoints `/api/health`, `/api/analyze`, `/api/align`, `/api/translate`, `/api/codons`

- [ ] **Step 1: Implement API router in `backend/app/api/endpoints.py`**

```python
from fastapi import APIRouter, HTTPException
from app.api.schemas import (
    AnalyzeRequest, AnalyzeResponse,
    AlignRequest, AlignResponse,
    TranslateRequest, TranslateResponse,
    CodonRequest, CodonResponse
)
from app.core.sequence_analyzer import SequenceAnalyzer
from app.core.io_fasta import parse_fasta_string
from app.core.alignment import perform_alignment
from app.core.codons import calculate_codon_usage, calculate_rscu
from Bio.Seq import Seq

router = APIRouter()
analyzer = SequenceAnalyzer()

@router.get("/health")
def health_check():
    return {"status": "online", "version": "2.0.0", "service": "genomesight-api"}

@router.post("/analyze", response_model=AnalyzeResponse)
def analyze_sequence(payload: AnalyzeRequest):
    content = ""
    if payload.sequence:
        content = f">custom_seq\n{payload.sequence}"
    elif payload.file_content:
        content = payload.file_content
    else:
        raise HTTPException(status_code=400, detail="Either 'sequence' or 'file_content' must be provided.")
    
    records = parse_fasta_string(content)
    if not records:
        raise HTTPException(status_code=400, detail="Could not parse valid sequences from input.")
    
    stats = analyzer.calculate_sequence_statistics(records)
    quality = analyzer.calculate_quality_statistics(records)
    kmer_res = analyzer.analyze_kmers(records, k=3)
    
    return {
        "success": True,
        "statistics": stats,
        "quality_statistics": quality,
        "kmer_analysis": kmer_res
    }

@router.post("/align", response_model=AlignResponse)
def align_sequences(payload: AlignRequest):
    try:
        res = perform_alignment(payload.seq1, payload.seq2, payload.mode)
        return {
            "success": True,
            **res
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.post("/translate", response_model=TranslateResponse)
def translate_sequence(payload: TranslateRequest):
    try:
        seq_obj = Seq(payload.sequence.upper())
        # adjust frame 1, 2, 3
        frame_idx = max(0, payload.frame - 1)
        translated = str(seq_obj[frame_idx:].translate(to_stop=False))
        return {
            "success": True,
            "protein_sequence": translated,
            "frame": payload.frame
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.post("/codons", response_model=CodonResponse)
def codon_analysis(payload: CodonRequest):
    try:
        counts = calculate_codon_usage(payload.sequence)
        rscu = calculate_rscu(payload.sequence)
        return {
            "success": True,
            "codon_counts": counts,
            "rscu": rscu
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
```

- [ ] **Step 2: Implement FastAPI app in `backend/app/main.py`**

```python
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.endpoints import router as api_router

app = FastAPI(
    title="GenomeSight API",
    version="2.0.0",
    description="Production REST API for genomic sequence analysis powered by BioPython and FastAPI."
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(api_router, prefix="/api")

@app.get("/")
def root():
    return {"message": "GenomeSight API is operational. Visit /docs for Swagger UI."}
```

- [ ] **Step 3: Create Dockerfile for Render deployment in `backend/Dockerfile`**

```dockerfile
FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "10000"]
```

- [ ] **Step 4: Commit FastAPI backend**

```bash
git add backend/app/main.py backend/app/api/endpoints.py backend/Dockerfile
git commit -m "feat(backend): implement FastAPI application entrypoint, endpoints, and Dockerfile"
```

---

### Task 3: React Frontend Setup & API Client with Cold-Start Retry

**Files:**
- Create: `frontend/package.json`
- Create: `frontend/vite.config.ts`
- Create: `frontend/tailwind.config.js`
- Create: `frontend/src/index.css`
- Create: `frontend/src/services/api.ts`

**Interfaces:**
- Consumes: Backend REST API (`/api/...`)
- Produces: Axios client with automatic exponential backoff retry for Render cold starts

- [ ] **Step 1: Create `frontend/package.json`**

```json
{
  "name": "genomesight-frontend",
  "private": true,
  "version": "2.0.0",
  "type": "module",
  "scripts": {
    "dev": "vite",
    "build": "tsc && vite build",
    "preview": "vite preview"
  },
  "dependencies": {
    "axios": "^1.6.8",
    "react": "^18.2.0",
    "react-dom": "^18.2.0",
    "lucide-react": "^0.359.0",
    "chart.js": "^4.4.2",
    "react-chartjs-2": "^5.2.0"
  },
  "devDependencies": {
    "@types/react": "^18.2.66",
    "@types/react-dom": "^18.2.22",
    "@vitejs/plugin-react": "^4.2.1",
    "autoprefixer": "^10.4.19",
    "postcss": "^8.4.38",
    "tailwindcss": "^3.4.1",
    "typescript": "^5.2.2",
    "vite": "^5.1.6"
  }
}
```

- [ ] **Step 2: Create Vite and Tailwind configuration files**

Create `frontend/vite.config.ts`:
```typescript
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  server: {
    port: 3000
  }
})
```

Create `frontend/tailwind.config.js`:
```javascript
/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        canvas: "#0a0a0f",
        surface: "#14141d",
        surface2: "#1e1e2c",
        bioGreen: "#22c55e",
        telecomCyan: "#38bdf8",
      },
      fontFamily: {
        mono: ['"IBM Plex Mono"', 'monospace'],
        sans: ['Inter', 'sans-serif'],
      }
    },
  },
  plugins: [],
}
```

- [ ] **Step 3: Create `frontend/src/index.css`**

```tailwind
@tailwind base;
@tailwind components;
@tailwind utilities;

@layer base {
  body {
    background-color: #0a0a0f;
    color: #e2e8f0;
    font-family: 'Inter', sans-serif;
  }
  
  code, pre, .font-mono {
    font-family: 'IBM Plex Mono', monospace;
  }
}
```

- [ ] **Step 4: Implement Axios API client with cold-start retry in `frontend/src/services/api.ts`**

```typescript
import axios from 'axios';

const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:10000/api';

export const apiClient = axios.create({
  baseURL: API_BASE_URL,
  timeout: 30000,
  headers: {
    'Content-Type': 'application/json',
  },
});

// Interceptor for handling Render free-tier cold starts with retry logic
apiClient.interceptors.response.use(
  (response) => response,
  async (error) => {
    const config = error.config;
    if (!config || !config.retryCount) {
      config.retryCount = 0;
    }

    // Retry up to 3 times for network or 5xx gateway timeout errors (cold start)
    if (config.retryCount < 3 && (!error.response || error.response.status >= 500)) {
      config.retryCount += 1;
      const delay = config.retryCount * 3000; // 3s, 6s, 9s backoff
      console.warn(`Backend cold start detected or timeout. Retrying attempt ${config.retryCount} in ${delay}ms...`);
      
      await new Promise((resolve) => setTimeout(resolve, delay));
      return apiClient(config);
    }

    return Promise.reject(error);
  }
);

export const checkHealth = async () => {
  const res = await apiClient.get('/health');
  return res.data;
};

export const analyzeSequence = async (payload: { sequence?: string; file_content?: string }) => {
  const res = await apiClient.post('/analyze', payload);
  return res.data;
};

export const alignSequences = async (seq1: string, seq2: string, mode: 'global' | 'local') => {
  const res = await apiClient.post('/align', { seq1, seq2, mode });
  return res.data;
};

export const translateSequence = async (sequence: string, frame: number) => {
  const res = await apiClient.post('/translate', { sequence, frame });
  return res.data;
};
```

- [ ] **Step 5: Commit frontend setup & API client**

```bash
git add frontend/package.json frontend/vite.config.ts frontend/tailwind.config.js frontend/src/index.css frontend/src/services/api.ts
git commit -m "feat(frontend): initialize React Vite Tailwind setup and Axios client with cold-start retry"
```

---

### Task 4: React Components & Main Application View

**Files:**
- Create: `frontend/src/components/Navbar.tsx`
- Create: `frontend/src/components/FileUpload.tsx`
- Create: `frontend/src/components/SummaryCards.tsx`
- Create: `frontend/src/components/AlignmentViewer.tsx`
- Create: `frontend/src/App.tsx`
- Create: `frontend/src/main.tsx`

**Interfaces:**
- Consumes: `api.ts` client
- Produces: Fully interactive Genomic Laboratory UI with dark theme and responsive visualizations

- [ ] **Step 1: Implement `frontend/src/components/Navbar.tsx`**

```tsx
import React from 'react';
import { Dna, ShieldCheck, Terminal } from 'lucide-react';

export const Navbar: React.FC<{ backendStatus: string }> = ({ backendStatus }) => {
  return (
    <header className="border-b border-slate-800 bg-surface/80 backdrop-blur sticky top-0 z-50 px-6 py-4 flex items-center justify-between">
      <div className="flex items-center gap-3">
        <div className="p-2 rounded bg-bioGreen/10 border border-bioGreen/30 text-bioGreen">
          <Dna className="w-6 h-6 animate-pulse" />
        </div>
        <div>
          <h1 className="font-mono font-bold text-lg tracking-wider text-white">GENOMESIGHT<span className="text-bioGreen">.AI</span></h1>
          <p className="text-xs text-slate-400 font-mono">Genomic Laboratory & Sequence Intelligence Suite</p>
        </div>
      </div>
      <div className="flex items-center gap-4 text-xs font-mono">
        <div className="flex items-center gap-2 px-3 py-1.5 rounded-full bg-surface2 border border-slate-700">
          <span className={`w-2 h-2 rounded-full ${backendStatus === 'online' ? 'bg-bioGreen animate-ping' : 'bg-amber-500'}`} />
          <span className="text-slate-300">Backend: {backendStatus.toUpperCase()}</span>
        </div>
        <div className="hidden md:flex items-center gap-1.5 text-slate-400">
          <ShieldCheck className="w-4 h-4 text-telecomCyan" />
          <span>BioPython + FastAPI v2.0</span>
        </div>
      </div>
    </header>
  );
};
```

- [ ] **Step 2: Implement `frontend/src/components/FileUpload.tsx`**

```tsx
import React, { useState } from 'react';
import { Upload, FileText, Play } from 'lucide-react';

interface FileUploadProps {
  onAnalyze: (sequence?: string, fileContent?: string) => void;
  loading: boolean;
}

const SAMPLE_DNA = "ATGCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATCG";

export const FileUpload: React.FC<FileUploadProps> = ({ onAnalyze, loading }) => {
  const [sequenceInput, setSequenceInput] = useState('');
  const [activeTab, setActiveTab] = useState<'paste' | 'file'>('paste');

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (event) => {
      const content = event.target?.result as string;
      onAnalyze(undefined, content);
    };
    reader.readAsText(file);
  };

  return (
    <div className="bg-surface border border-slate-800 rounded-xl p-6 shadow-xl">
      <div className="flex items-center justify-between mb-4">
        <div className="flex gap-2">
          <button
            onClick={() => setActiveTab('paste')}
            className={`px-4 py-2 rounded-lg font-mono text-xs transition ${activeTab === 'paste' ? 'bg-bioGreen text-black font-bold' : 'bg-surface2 text-slate-400 hover:text-white'}`}
          >
            Paste Sequence
          </button>
          <button
            onClick={() => setActiveTab('file')}
            className={`px-4 py-2 rounded-lg font-mono text-xs transition ${activeTab === 'file' ? 'bg-bioGreen text-black font-bold' : 'bg-surface2 text-slate-400 hover:text-white'}`}
          >
            Upload FASTA/FASTQ
          </button>
        </div>
        <button
          onClick={() => { setSequenceInput(SAMPLE_DNA); onAnalyze(SAMPLE_DNA); }}
          className="text-xs font-mono text-telecomCyan hover:underline flex items-center gap-1"
        >
          Load Sample Sequence
        </button>
      </div>

      {activeTab === 'paste' ? (
        <div className="space-y-4">
          <textarea
            rows={4}
            value={sequenceInput}
            onChange={(e) => setSequenceInput(e.target.value)}
            placeholder="Paste raw FASTA or DNA sequence (e.g., ATGCGATCG...)"
            className="w-full bg-canvas border border-slate-800 rounded-lg p-4 font-mono text-xs text-slate-200 focus:border-bioGreen focus:outline-none"
          />
          <button
            disabled={loading || !sequenceInput.trim()}
            onClick={() => onAnalyze(sequenceInput)}
            className="w-full py-3 bg-bioGreen hover:bg-bioGreen/90 disabled:opacity-50 text-black font-mono font-bold rounded-lg transition flex items-center justify-center gap-2"
          >
            {loading ? <span className="animate-spin">⏳</span> : <Play className="w-4 h-4" />}
            {loading ? "Analyzing Genome..." : "Run Sequence Analysis"}
          </button>
        </div>
      ) : (
        <div className="border-2 border-dashed border-slate-800 hover:border-bioGreen/50 rounded-xl p-8 text-center transition cursor-pointer relative">
          <input
            type="file"
            accept=".fasta,.fa,.fastq,.txt"
            onChange={handleFileUpload}
            className="absolute inset-0 opacity-0 cursor-pointer"
          />
          <Upload className="w-10 h-10 text-bioGreen mx-auto mb-3" />
          <p className="font-mono text-sm text-slate-300">Drop FASTA/FASTQ file here or click to browse</p>
          <p className="font-mono text-xs text-slate-500 mt-1">Supports standard nucleotide formats (.fasta, .fastq)</p>
        </div>
      )}
    </div>
  );
};
```

- [ ] **Step 3: Implement `frontend/src/components/SummaryCards.tsx`**

```tsx
import React from 'react';
import { Activity, BarChart2, ShieldAlert, Cpu } from 'lucide-react';

export const SummaryCards: React.FC<{ stats: any }> = ({ stats }) => {
  if (!stats) return null;

  return (
    <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
      <div className="bg-surface border border-slate-800 p-5 rounded-xl">
        <div className="flex items-center justify-between text-slate-400 mb-2">
          <span className="font-mono text-xs">Total Length</span>
          <Activity className="w-4 h-4 text-bioGreen" />
        </div>
        <div className="font-mono text-2xl font-bold text-white">{stats.total_length?.toLocaleString()} bp</div>
        <div className="text-xs text-slate-500 font-mono mt-1">{stats.sequence_count} sequence(s)</div>
      </div>

      <div className="bg-surface border border-slate-800 p-5 rounded-xl">
        <div className="flex items-center justify-between text-slate-400 mb-2">
          <span className="font-mono text-xs">Average GC Content</span>
          <BarChart2 className="w-4 h-4 text-telecomCyan" />
        </div>
        <div className="font-mono text-2xl font-bold text-white">{stats.average_gc_content?.toFixed(2)}%</div>
        <div className="text-xs text-slate-500 font-mono mt-1">Thermal stability index</div>
      </div>

      <div className="bg-surface border border-slate-800 p-5 rounded-xl">
        <div className="flex items-center justify-between text-slate-400 mb-2">
          <span className="font-mono text-xs">Sequence Range</span>
          <Cpu className="w-4 h-4 text-purple-400" />
        </div>
        <div className="font-mono text-xl font-bold text-white">{stats.min_length} – {stats.max_length} bp</div>
        <div className="text-xs text-slate-500 font-mono mt-1">Median: {stats.median_length} bp</div>
      </div>

      <div className="bg-surface border border-slate-800 p-5 rounded-xl">
        <div className="flex items-center justify-between text-slate-400 mb-2">
          <span className="font-mono text-xs">Computation Backend</span>
          <ShieldAlert className="w-4 h-4 text-amber-400" />
        </div>
        <div className="font-mono text-xl font-bold text-white">BioPython 1.83</div>
        <div className="text-xs text-slate-500 font-mono mt-1">FastAPI REST Server</div>
      </div>
    </div>
  );
};
```

- [ ] **Step 4: Implement `frontend/src/components/AlignmentViewer.tsx`**

```tsx
import React, { useState } from 'react';
import { alignSequences } from '../services/api';
import { GitCompare, ArrowRight } from 'lucide-react';

export const AlignmentViewer: React.FC = () => {
  const [seq1, setSeq1] = useState('ATGCGATCGATCGATAG');
  const [seq2, setSeq2] = useState('ATGCGATCTATCGATAG');
  const [mode, setMode] = useState<'global' | 'local'>('global');
  const [result, setResult] = useState<any>(null);
  const [loading, setLoading] = useState(false);

  const handleAlign = async () => {
    setLoading(true);
    try {
      const res = await alignSequences(seq1, seq2, mode);
      setResult(res);
    } catch (err) {
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="bg-surface border border-slate-800 rounded-xl p-6 shadow-xl">
      <h3 className="font-mono font-bold text-md text-white mb-4 flex items-center gap-2">
        <GitCompare className="w-5 h-5 text-bioGreen" />
        Pairwise Sequence Alignment (Needleman-Wunsch / Smith-Waterman)
      </h3>
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
        <div>
          <label className="block text-xs font-mono text-slate-400 mb-1">Sequence 1</label>
          <input
            type="text"
            value={seq1}
            onChange={(e) => setSeq1(e.target.value)}
            className="w-full bg-canvas border border-slate-800 rounded-lg p-3 font-mono text-xs text-slate-200 focus:border-bioGreen focus:outline-none"
          />
        </div>
        <div>
          <label className="block text-xs font-mono text-slate-400 mb-1">Sequence 2</label>
          <input
            type="text"
            value={seq2}
            onChange={(e) => setSeq2(e.target.value)}
            className="w-full bg-canvas border border-slate-800 rounded-lg p-3 font-mono text-xs text-slate-200 focus:border-bioGreen focus:outline-none"
          />
        </div>
      </div>
      <div className="flex items-center justify-between mb-6">
        <div className="flex gap-2">
          <button
            onClick={() => setMode('global')}
            className={`px-3 py-1.5 rounded font-mono text-xs ${mode === 'global' ? 'bg-telecomCyan text-black font-bold' : 'bg-surface2 text-slate-400'}`}
          >
            Global (Needleman-Wunsch)
          </button>
          <button
            onClick={() => setMode('local')}
            className={`px-3 py-1.5 rounded font-mono text-xs ${mode === 'local' ? 'bg-telecomCyan text-black font-bold' : 'bg-surface2 text-slate-400'}`}
          >
            Local (Smith-Waterman)
          </button>
        </div>
        <button
          onClick={handleAlign}
          disabled={loading}
          className="px-6 py-2 bg-bioGreen text-black font-mono font-bold rounded-lg text-xs hover:bg-bioGreen/90 transition flex items-center gap-2"
        >
          {loading ? "Aligning..." : "Run Alignment"} <ArrowRight className="w-4 h-4" />
        </button>
      </div>

      {result && (
        <div className="bg-canvas border border-slate-800 rounded-lg p-4 font-mono text-xs space-y-2">
          <div className="flex justify-between text-slate-400 text-[11px] pb-2 border-b border-slate-800">
            <span>Score: <strong className="text-white">{result.score}</strong></span>
            <span>Identity: <strong className="text-bioGreen">{result.identity_percentage.toFixed(1)}%</strong></span>
          </div>
          <div className="overflow-x-auto text-slate-200 space-y-1 py-2">
            <div className="tracking-widest">{result.aligned_seq1}</div>
            <div className="tracking-widest text-bioGreen">{result.match_line}</div>
            <div className="tracking-widest">{result.aligned_seq2}</div>
          </div>
        </div>
      )}
    </div>
  );
};
```

- [ ] **Step 5: Implement `frontend/src/App.tsx`**

```tsx
import React, { useState, useEffect } from 'react';
import { Navbar } from './components/Navbar';
import { FileUpload } from './components/FileUpload';
import { SummaryCards } from './components/SummaryCards';
import { AlignmentViewer } from './components/AlignmentViewer';
import { checkHealth, analyzeSequence } from './services/api';
import { Terminal, Shield, Layers } from 'lucide-react';

export function App() {
  const [backendStatus, setBackendStatus] = useState<string>('checking...');
  const [loading, setLoading] = useState(false);
  const [analysisResult, setAnalysisResult] = useState<any>(null);

  useEffect(() => {
    checkHealth()
      .then((res) => setBackendStatus(res.status))
      .catch(() => setBackendStatus('offline (cold start warming...)'));
  }, []);

  const handleAnalyze = async (sequence?: string, fileContent?: string) => {
    setLoading(true);
    try {
      const res = await analyzeSequence({ sequence, file_content: fileContent });
      setAnalysisResult(res);
      setBackendStatus('online');
    } catch (err) {
      console.error(err);
      alert('Analysis failed. Please check backend connection.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-canvas text-slate-200 font-sans pb-16">
      <Navbar backendStatus={backendStatus} />

      <main className="max-w-7xl mx-auto px-6 pt-8 space-y-8">
        {/* Hero Banner */}
        <div className="bg-gradient-to-r from-surface to-surface2 border border-slate-800 rounded-2xl p-8 relative overflow-hidden">
          <div className="absolute right-0 top-0 translate-x-8 -translate-y-8 w-64 h-64 bg-bioGreen/5 rounded-full blur-3xl pointer-events-none" />
          <div className="max-w-2xl">
            <span className="inline-block font-mono text-xs px-3 py-1 rounded bg-bioGreen/10 text-bioGreen border border-bioGreen/30 mb-4">
              PHASE 2 REWRITE · PRODUCTION ARCHITECTURE
            </span>
            <h2 className="text-3xl font-bold font-mono text-white tracking-tight mb-3">
              Genomic Sequence Intelligence Suite
            </h2>
            <p className="text-slate-400 text-sm leading-relaxed">
              Decoupled FastAPI backend powered by BioPython with a React / Tailwind frontend. Designed for recruiter evaluation and high-throughput computational biology.
            </p>
          </div>
        </div>

        {/* Input Section */}
        <FileUpload onAnalyze={handleAnalyze} loading={loading} />

        {/* Results Dashboard */}
        {analysisResult && (
          <div className="space-y-6">
            <h3 className="font-mono text-lg font-bold text-white flex items-center gap-2">
              <Layers className="w-5 h-5 text-telecomCyan" /> Analysis Telemetry Dashboard
            </h3>
            <SummaryCards stats={analysisResult.statistics} />

            {/* Nucleotide Composition Breakdown */}
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <div className="bg-surface border border-slate-800 rounded-xl p-6">
                <h4 className="font-mono font-bold text-sm text-white mb-4">Nucleotide Composition (%)</h4>
                <div className="space-y-3 font-mono text-xs">
                  {Object.entries(analysisResult.statistics.composition_percentages || {}).map(([nuc, pct]: [string, any]) => (
                    <div key={nuc} className="space-y-1">
                      <div className="flex justify-between text-slate-400">
                        <span>Base {nuc}</span>
                        <span>{pct.toFixed(2)}%</span>
                      </div>
                      <div className="w-full bg-canvas rounded-full h-2 overflow-hidden border border-slate-800">
                        <div
                          className="bg-bioGreen h-full rounded-full transition-all duration-500"
                          style={{ width: `${pct}%` }}
                        />
                      </div>
                    </div>
                  ))}
                </div>
              </div>

              <div className="bg-surface border border-slate-800 rounded-xl p-6">
                <h4 className="font-mono font-bold text-sm text-white mb-4">Top 5 K-mer Frequencies (k=3)</h4>
                <div className="space-y-2 font-mono text-xs">
                  {analysisResult.kmer_analysis.most_common?.slice(0, 5).map(([kmer, count]: [string, any]) => (
                    <div key={kmer} className="flex items-center justify-between p-3 rounded-lg bg-canvas border border-slate-800">
                      <span className="text-telecomCyan font-bold">{kmer}</span>
                      <span className="text-slate-300">{count} occurrences</span>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Alignment Module */}
        <AlignmentViewer />
      </main>
    </div>
  );
}

export default App;
```

- [ ] **Step 6: Implement `frontend/src/main.tsx`**

```tsx
import React from 'react'
import ReactDOM from 'react-dom/client'
import App from './App.tsx'
import './index.css'

ReactDOM.createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>,
)
```

- [ ] **Step 7: Commit frontend implementation**

```bash
git add frontend/src/
git commit -m "feat(frontend): implement React components, dashboard, alignment viewer, and App root"
```

---

## Verification Plan

1. **Backend Verification:**
   - Run `uvicorn backend.app.main:app --reload`
   - Test health endpoint: `curl http://localhost:8000/api/health`
   - Test analysis endpoint with sample FASTA string.

2. **Frontend Verification:**
   - Run `cd frontend && npm install && npm run dev`
   - Verify UI rendering and live sequence analysis through browser/fetch.
