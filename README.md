![Project Banner](assets/banner.png)
# GenomeSight

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

GenomeSight is an interactive genome sequence analysis toolkit implemented in Python using Streamlit. It supports exploratory analysis of DNA and RNA sequences, including GC content calculation, k-mer profiling, ORF detection, motif search, restriction site identification, and quality assessment with interactive visualizations.

🌐 Live Demo: https://genomesight.arditmishra.com

## Features

- **Sequence Statistics**: Calculate GC content, nucleotide composition, and sequence complexity
- **K-mer Analysis**: Frequency analysis of subsequences (configurable k-mer size), with an optional native (Cython) counting accelerator for large genomes — see [Native k-mer accelerator](#native-k-mer-accelerator)
- **ORF Detection**: Find Open Reading Frames with customizable minimum length
- **Codon Usage**: Per-ORF and whole-sequence codon counts with relative synonymous codon usage (RSCU)
- **Motif Search**: Pattern matching with IUPAC ambiguity code support
- **Restriction Sites**: Search for common restriction enzyme recognition sequences
- **Quality Analysis**: FASTQ quality score visualization and statistics
- **Interactive Plots**: Dynamic Plotly visualizations with dark theme
- **Export Options**: Download results as JSON, CSV, or comprehensive ZIP archives

## Installation

### Prerequisites

- Python 3.11 or higher
- pip package manager

### Quick Start

```bash
# Clone the repository
git clone https://github.com/Ardit-Mishra/genomesight.git
cd genomesight

# (Optional) Create a virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies (recommended)
pip install -e .

# Run the application
streamlit run main.py --server.port 5000
```

The application will be available at `http://localhost:5000`

## Usage

### Basic Analysis

1. **Upload a sequence file** (FASTA, FASTQ, or GenBank format)
2. View automatic analysis results including:
   - Sequence count and total length
   - Average GC content
   - Nucleotide composition chart
   - Length distribution

### Advanced Analysis

Use the sidebar to access additional tools:

- **K-mer Analysis**: Set k-mer size (2-10) and analyze frequency patterns
- **ORF Detection**: Find protein-coding regions with adjustable minimum length
- **Codon Usage**: Codon counts and RSCU, both per detected ORF (aggregated) and for the raw uploaded sequence
- **Motif Search**: Enter IUPAC patterns or select restriction enzymes

### Supported File Formats

| Format | Extensions | Description |
|--------|------------|-------------|
| FASTA | `.fasta`, `.fa`, `.fna` | Standard sequence format |
| FASTQ | `.fastq`, `.fq` | Sequences with quality scores |
| GenBank | `.gb`, `.gbk`, `.genbank` | Sequence plus structured feature annotations (source, CDS, gene, etc.) |

### Example

```python
# Using the core library programmatically
from app.core.sequence_analyzer import SequenceAnalyzer
from app.core.io_fasta import parse_fasta

# Parse sequences
with open('sample_data/example.fasta', 'r') as f:
    sequences = parse_fasta(f.read())

# Analyze
analyzer = SequenceAnalyzer()
for seq in sequences:
    gc = analyzer.calculate_gc_content(str(seq.seq))
    print(f"{seq.id}: {gc:.2f}% GC content")
```

## Project Structure

```
genomesight/
├── main.py                     # Streamlit application entry point
├── app/
│   ├── core/                   # Core analysis modules (standalone library)
│   │   ├── sequence_analyzer.py    # GC content, composition, k-mers
│   │   ├── kmer_native.py          # K-mer counting dispatcher (native/Python)
│   │   ├── _kmer_accel.pyx         # Optional Cython k-mer accelerator source
│   │   ├── io_fasta.py             # FASTA/FASTQ parsing
│   │   ├── validation.py           # Sequence validation
│   │   ├── orf.py                  # ORF detection
│   │   ├── motifs.py               # Pattern matching
│   │   ├── alignment.py            # Pairwise alignment
│   │   ├── plots.py                # Plotly visualizations
│   │   └── export.py               # Report generation
│   └── ui/                     # Streamlit UI components
│       ├── components.py           # Reusable widgets
│       └── styles.py               # CSS theming
├── benchmarks/                     # Native vs Python k-mer timing script + results
├── build_kmer_ext.py               # Optional: compiles the Cython k-mer accelerator
├── .streamlit/                     # Streamlit configuration
├── docs/                           # Documentation and GitHub instructions
├── tests/                          # Unit tests
├── sample_data/                    # Example FASTA/FASTQ files
├── pyproject.toml                  # Project configuration and dependencies
└── uv.lock                         # Locked dependency versions
```

## API Reference

### SequenceAnalyzer

```python
from app.core.sequence_analyzer import SequenceAnalyzer

analyzer = SequenceAnalyzer()

# GC Content (returns percentage 0-100)
gc = analyzer.calculate_gc_content("ATGCATGC")  # 50.0

# Nucleotide Composition
comp = analyzer.calculate_nucleotide_composition("ATGC")
# {'A': 1, 'T': 1, 'C': 1, 'G': 1, 'N': 0}

# Reverse Complement
rc = analyzer.reverse_complement("ATGC")  # "GCAT"

# K-mer Analysis -- uses the native accelerator when it's built, otherwise
# an equivalent pure-Python counter; either way the result includes which
# backend actually ran.
kmers = analyzer.analyze_kmers(sequences, k=3)
kmers['backend']  # "native" or "python"
```

### ORF Detection

```python
from app.core.orf import find_orfs

orfs = find_orfs(
    sequence="ATGAAAAAAAAATGA",
    min_length=9,
    include_reverse=True
)

for orf in orfs:
    print(f"ORF: {orf.start}-{orf.end}, Length: {orf.length_nt} bp")
```

### Motif Search

```python
from app.core.motifs import search_motif, RESTRICTION_ENZYMES

# Search with IUPAC pattern
matches = search_motif(sequences, pattern="ATGR")  # R = A or G

# Search for restriction site
matches = search_motif(sequences, pattern=RESTRICTION_ENZYMES['EcoRI'])
```

## Native k-mer accelerator

K-mer counting (the sliding-window scan over every sequence for a chosen
`k`) has two implementations:

1. **Pure Python** (`app/core/kmer_native._count_kmers_python`) — always
   available, zero extra dependencies, what the app uses out of the box.
2. **Native, Cython-compiled** (`app/core/_kmer_accel.pyx`) — an optional
   compiled extension for faster counting on large genomes. It is **not**
   part of the app's required dependencies and is **not** committed as a
   binary; you build it yourself if you want it.

`app/core/kmer_native.count_kmers()` tries the compiled extension first
and transparently falls back to the Python counter if it isn't present —
the app never crashes for lack of a compiler. The two implementations are
required to return byte-for-byte identical counts
(`tests/test_kmer_native.py` checks this on hand-picked and randomized
inputs whenever the extension is built), and the Streamlit UI's K-mer
Analysis panel always states outright which backend produced the numbers
you're looking at — it is never left unstated.

### Building it

Requires a C compiler on the host (MSVC on Windows, gcc/clang on
Linux/macOS) plus Cython and setuptools — build-time only, not runtime
dependencies:

```bash
uv run --python 3.12 --with cython --with setuptools \
    python build_kmer_ext.py build_ext --inplace
```

or, with Cython and a compiler already on `PATH`:

```bash
pip install -e ".[native]"
python build_kmer_ext.py build_ext --inplace
```

On success this drops a compiled module next to `_kmer_accel.pyx` (e.g.
`_kmer_accel.cp312-win_amd64.pyd` on Windows or
`_kmer_accel.cpython-312-x86_64-linux-gnu.so` on Linux); `kmer_native.py`
picks it up automatically on the next import, no other wiring needed. If
the build fails or is skipped — no compiler on the deploy host, for
instance — the app still runs correctly on the Python path and says so in
the UI.

### Benchmarking it

```bash
uv run --python 3.12 python benchmarks/benchmark_kmer.py
```

Times both paths (when native is built) across several sequence lengths
and `k` values on identical input, and writes
`benchmarks/kmer_benchmark_results.json`. It refuses to report a speedup
if the two paths ever disagree on the result.

Measured on this machine (Windows, Python 3.12, MSVC 14.50 build, min of
3 runs per cell, random sequence):

| sequence length | k | Python | Native | speedup |
|---|---|---|---|---|
| 1,000 | 4 | 0.86 ms | 0.12 ms | 7.0x |
| 1,000 | 8 | 1.07 ms | 0.12 ms | 9.3x |
| 1,000 | 12 | 1.22 ms | 0.12 ms | 10.3x |
| 10,000 | 4 | 7.55 ms | 1.25 ms | 6.0x |
| 10,000 | 8 | 10.93 ms | 1.20 ms | 9.1x |
| 10,000 | 12 | 12.78 ms | 1.17 ms | 11.0x |
| 100,000 | 4 | 76.70 ms | 12.60 ms | 6.1x |
| 100,000 | 8 | 131.43 ms | 16.98 ms | 7.7x |
| 100,000 | 12 | 307.36 ms | 51.72 ms | 5.9x |
| 500,000 | 4 | 634.01 ms | 103.88 ms | 6.1x |
| 500,000 | 8 | 1,180.37 ms | 303.17 ms | 3.9x |
| 500,000 | 12 | 1,323.71 ms | 343.17 ms | 3.9x |

**The speedup shrinks as the input grows**, from ~11x at 10 kb down to
~3.9x at 500 kb. That is the honest and expected shape: the Cython loop
removes Python's per-character interpreter overhead, but both paths still
build the same Python dictionary, and at large `k` on a long sequence
nearly every window is a distinct key — so dictionary growth, hashing and
memory dominate the runtime and neither implementation can avoid them.
The accelerator makes the *scan* cheap; it does not make the *counting*
cheap. Anyone quoting a single headline multiple for this is quoting the
best cell in the table.

The exact numbers above are regenerated by the command in this section
and live in `benchmarks/kmer_benchmark_results.json`, whose `environment`
block records whether the native path was actually available for the run
— so a table produced without a compiler is self-evidently a Python-only
table rather than a silent half-measurement.

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_sequence_analyzer.py -v

# Run with coverage
pytest tests/ --cov=app --cov-report=html
```
## Core Dependencies

| Package | Purpose |
|---------|---------|
| Biopython | Sequence analysis and file parsing |
| NumPy | Numerical computation |
| Pandas | Data handling |
| Plotly | Interactive visualization |
| Streamlit | Web application framework |


## Citation

If you use this software in your research, please cite:

```bibtex
@software{mishra2025genomesight,
  author = {Mishra, Ardit},
  title = {GenomeSight: Interactive Genome Sequence Analysis Toolkit},
  year = {2025},
  url = {https://github.com/Ardit-Mishra/genomesight},
  version = {1.0.0}
}
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Ardit Mishra**

- GitHub: [@Ardit-Mishra](https://github.com/Ardit-Mishra)
- Website: https://arditmishra.com
  
## Acknowledgments

- [Biopython](https://biopython.org/) for sequence analysis tools
- [Streamlit](https://streamlit.io/) for the web framework
- [Plotly](https://plotly.com/) for interactive visualizations
