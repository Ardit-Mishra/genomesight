"""
Core analysis modules for the Genome Sequencing Analyzer.

This package contains all pure Python functions and classes for sequence analysis.
These modules do not depend on Streamlit and can be used independently.

Modules:
    - sequence_analyzer: Core sequence analysis functions
    - io_fasta: FASTA/FASTQ file parsing and export
    - validation: Sequence validation and type detection
    - orf: Open Reading Frame detection
    - motifs: Motif discovery and IUPAC pattern matching
    - alignment: Pairwise sequence alignment
    - plots: Plotly visualization generators
    - export: Report generation and export utilities
"""

from .sequence_analyzer import SequenceAnalyzer
from .io_fasta import parse_fasta, parse_fastq, export_to_fasta
from .validation import validate_sequence, detect_sequence_type
from .orf import find_orfs, translate_sequence
from .motifs import search_motif, IUPAC_CODES
from .plots import (
    create_gc_plot,
    create_composition_plot,
    create_kmer_plot,
    create_quality_plot
)
