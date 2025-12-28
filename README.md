# Genome Sequencing Analyzer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.0+-FF4B4B.svg)](https://streamlit.io/)
[![Biopython](https://img.shields.io/badge/Biopython-1.80+-green.svg)](https://biopython.org/)

A comprehensive bioinformatics web application for DNA/RNA sequence analysis. Built with Python and Streamlit, featuring GC content analysis, nucleotide composition, k-mer frequency analysis, Open Reading Frame (ORF) detection, motif discovery, and interactive visualizations.

<!-- Screenshot placeholder - add screenshot after deployment
<p align="center">
  <img src="assets/screenshot.png" alt="Genome Sequencing Analyzer Screenshot" width="800">
</p>
-->

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
git clone https://github.com/yourusername/genome-sequencing-analyzer.git
cd genome-sequencing-analyzer

# Install dependencies
pip install -r requirements.txt

# Run the application
streamlit run main.py --server.port 5000
```

The application will be available at `http://localhost:5000`

### Using Conda

```bash
# Create environment
conda create -n genome-analyzer python=3.11
conda activate genome-analyzer

# Install dependencies
pip install -r requirements.txt

# Run
streamlit run main.py --server.port 5000
```

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
genome-sequencing-analyzer/
├── main.py                     # Application entry point
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
├── tests/                      # Unit tests
├── sample_data/                # Example files
└── requirements.txt            # Dependencies
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

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| Biopython | >=1.80 | Sequence analysis and file parsing |
| NumPy | >=1.24 | Numerical computations |
| Pandas | >=2.0 | Data manipulation |
| Plotly | >=5.15 | Interactive visualizations |
| Streamlit | >=1.28 | Web application framework |

## Citation

If you use this software in your research, please cite:

```bibtex
@software{mishra2025genome,
  author = {Mishra, Ardit},
  title = {Genome Sequencing Analyzer: A Comprehensive DNA/RNA Analysis Tool},
  year = {2025},
  url = {https://github.com/yourusername/genome-sequencing-analyzer},
  version = {1.0.0}
}
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Ardit Mishra**

- GitHub: [@yourusername](https://github.com/yourusername)

## Acknowledgments

- [Biopython](https://biopython.org/) for sequence analysis tools
- [Streamlit](https://streamlit.io/) for the web framework
- [Plotly](https://plotly.com/) for interactive visualizations
