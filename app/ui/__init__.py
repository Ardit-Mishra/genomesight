"""
UI Components for the Genome Sequencing Analyzer.

This package contains reusable Streamlit UI components and styling utilities.

Modules:
    - components: Reusable UI widgets and display functions
    - styles: CSS styling and theme configuration
"""

from .components import (
    display_sequence_card,
    display_metrics_row,
    display_analysis_results,
    file_uploader_section
)
from .styles import apply_dark_theme, get_nucleotide_colors
