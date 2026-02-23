# GenomeSight

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

GenomeSight is an interactive genome sequence analysis toolkit implemented in Python using Streamlit. It supports exploratory analysis of DNA and RNA sequences, including GC content calculation, k-mer profiling, ORF detection, motif search, restriction site identification, and quality assessment with interactive visualizations.

🌐 Live Demo: https://genomesight.arditmishra.com

## Features

- **Sequence Statistics**: Calculate GC content, nucleotide composition, and sequence complexity
- **K-mer Analysis**: Frequency analysis of subsequences (configurable k-mer size)
- **ORF Detection**: Find Open Reading Frames with customizable minimum length
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

1. **Upload a sequence file** (FASTA or FASTQ format)
2. View automatic analysis results including:
   - Sequence count and total length
   - Average GC content
   - Nucleotide composition chart
   - Length distribution

### Advanced Analysis

Use the sidebar to access additional tools:

- **K-mer Analysis**: Set k-mer size (2-10) and analyze frequency patterns
- **ORF Detection**: Find protein-coding regions with adjustable minimum length
- **Motif Search**: Enter IUPAC patterns or select restriction enzymes

### Supported File Formats

| Format | Extensions | Description |
|--------|------------|-------------|
| FASTA | `.fasta`, `.fa`, `.fna` | Standard sequence format |
| FASTQ | `.fastq`, `.fq` | Sequences with quality scores |

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

# K-mer Analysis
kmers = analyzer.analyze_kmers(sequences, k=3)
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
