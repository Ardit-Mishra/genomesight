# Genome Sequencing Analyzer

A professional bioinformatics web application for comprehensive DNA/RNA sequence analysis. Built with Streamlit, featuring GC content analysis, nucleotide composition, k-mer frequency, ORF detection, motif discovery, and interactive visualizations.

## Quick Start

```bash
streamlit run main.py --server.port 5000
```

## Project Structure

```
genome-sequencing-analyzer/
├── main.py                     # Application entry point
├── app/
│   ├── __init__.py
│   ├── core/                   # Core analysis modules (no Streamlit deps)
│   │   ├── __init__.py
│   │   ├── sequence_analyzer.py    # Main analysis class
│   │   ├── io_fasta.py             # FASTA/FASTQ parsing
│   │   ├── validation.py           # Sequence validation
│   │   ├── orf.py                  # ORF detection
│   │   ├── motifs.py               # Motif discovery
│   │   ├── alignment.py            # Pairwise alignment
│   │   ├── plots.py                # Plotly visualizations
│   │   └── export.py               # Report generation
│   └── ui/                     # Streamlit UI components
│       ├── __init__.py
│       ├── components.py           # Reusable widgets
│       └── styles.py               # CSS and theming
├── tests/                      # Unit tests
│   ├── __init__.py
│   ├── test_sequence_analyzer.py
│   ├── test_io_fasta.py
│   └── test_orf.py
├── assets/                     # Static assets
└── .streamlit/
    └── config.toml             # Streamlit configuration
```

## Architecture

### Core Layer (`app/core/`)
Pure Python modules with no Streamlit dependencies. Can be used as a standalone library.

- **sequence_analyzer.py**: GC content, composition, complexity, k-mer analysis
- **io_fasta.py**: FASTA/FASTQ parsing and export
- **validation.py**: Sequence validation with IUPAC support
- **orf.py**: Open Reading Frame detection and translation
- **motifs.py**: Pattern matching with IUPAC codes, restriction sites
- **alignment.py**: Pairwise sequence alignment
- **plots.py**: Plotly chart generators with dark theme
- **export.py**: JSON, CSV, and report generation

### UI Layer (`app/ui/`)
Streamlit-specific components and styling.

- **components.py**: Reusable display widgets
- **styles.py**: CSS styling and color palette

## Design System

### Color Palette
- Background: #1A1D29 (dark blue-gray)
- Surface: #25293C (charcoal)
- Primary: #00C9A7 (teal)
- Text: #E8EAED (light gray)
- DNA bases: A=#FF6B9D, T=#FFC252, C=#00C9A7, G=#4D96FF

## Dependencies

- **Biopython**: Sequence analysis and file parsing
- **NumPy**: Numerical computations
- **Pandas**: Data manipulation
- **Plotly**: Interactive visualizations
- **Streamlit**: Web framework

## Testing

```bash
pytest tests/ -v
```

## Author

Ardit Mishra

## Version

1.0.0
