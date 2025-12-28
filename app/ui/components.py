"""
Reusable UI Components Module

This module provides reusable Streamlit UI components for displaying
sequence data, analysis results, and interactive elements.

Components:
    - display_sequence_card: Formatted sequence display
    - display_metrics_row: Row of metric cards
    - display_analysis_results: Complete results panel
    - file_uploader_section: Styled file upload area
"""

import streamlit as st
from typing import Dict, Any, List, Optional
from .styles import COLORS, NUCLEOTIDE_COLORS, format_sequence_html


def display_sequence_card(
    sequence: str,
    title: str = "Sequence",
    max_display: int = 200
) -> None:
    """
    Display a sequence in a formatted card.
    
    Args:
        sequence: Nucleotide sequence string
        title: Card title
        max_display: Maximum bases to display
    """
    display_seq = sequence[:max_display]
    if len(sequence) > max_display:
        display_seq += f"... ({len(sequence) - max_display} more bases)"
    
    st.markdown(f"""
    <div class="metric-card">
        <h4 style="color: {COLORS['primary']}; margin-bottom: 0.5rem;">{title}</h4>
        <p style="font-family: monospace; word-break: break-all; color: {COLORS['text_primary']};">
            {format_sequence_html(display_seq)}
        </p>
        <p style="color: {COLORS['text_secondary']}; font-size: 0.9rem;">
            Length: {len(sequence):,} bp
        </p>
    </div>
    """, unsafe_allow_html=True)


def display_metrics_row(metrics: Dict[str, Any]) -> None:
    """
    Display a row of metrics.
    
    Args:
        metrics: Dictionary of metric name to value
    """
    cols = st.columns(len(metrics))
    for col, (name, value) in zip(cols, metrics.items()):
        with col:
            if isinstance(value, float):
                col.metric(name, f"{value:.2f}")
            else:
                col.metric(name, value)


def display_analysis_results(results: Dict[str, Any]) -> None:
    """
    Display comprehensive analysis results.
    
    Args:
        results: Analysis results dictionary
    """
    if not results:
        st.info("No results to display")
        return
    
    stats = results.get('basic_stats', {})
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Sequences", stats.get('sequence_count', 0))
    with col2:
        st.metric("Total Length", f"{stats.get('total_length', 0):,} bp")
    with col3:
        st.metric("Avg GC%", f"{stats.get('average_gc_content', 0):.1f}%")
    with col4:
        st.metric("Avg Length", f"{stats.get('average_length', 0):.0f} bp")


def file_uploader_section() -> None:
    """
    Display the styled file upload section header.
    """
    st.markdown(f"""
    <div style="background: linear-gradient(135deg, {COLORS['surface']} 0%, #1E2330 100%); 
                border: 2px dashed {COLORS['primary']}; 
                border-radius: 12px; 
                padding: 2.5rem; 
                text-align: center; 
                margin: 1.5rem 0;
                box-shadow: 0 4px 6px rgba(0, 201, 167, 0.1);">
        <div style="font-size: 3rem; margin-bottom: 1rem;">🧬</div>
        <h3 style="color: {COLORS['primary']}; margin: 0 0 0.5rem 0; font-size: 1.5rem;">
            Upload Your Sequence Files
        </h3>
        <p style="color: {COLORS['text_secondary']}; margin: 0 0 1rem 0;">
            Drag and drop files or click to browse
        </p>
        <p style="color: {COLORS['text_primary']}; font-size: 0.9rem;">
            Supported: FASTA (.fasta, .fa) • FASTQ (.fastq, .fq)
        </p>
    </div>
    """, unsafe_allow_html=True)


def display_orf_table(orfs: List) -> None:
    """
    Display ORF results in a formatted table.
    
    Args:
        orfs: List of ORF objects
    """
    if not orfs:
        st.info("No ORFs found")
        return
    
    import pandas as pd
    
    data = []
    for i, orf in enumerate(orfs[:50], 1):
        data.append({
            '#': i,
            'Frame': f"{orf.strand}{orf.frame}",
            'Start': orf.start,
            'End': orf.end,
            'Length (nt)': orf.length_nt,
            'Length (aa)': orf.length_aa,
            'GC%': f"{orf.gc_content:.1f}"
        })
    
    df = pd.DataFrame(data)
    st.dataframe(df, use_container_width=True)


def display_motif_results(matches: List) -> None:
    """
    Display motif search results.
    
    Args:
        matches: List of MotifMatch objects
    """
    if not matches:
        st.info("No matches found")
        return
    
    st.success(f"Found {len(matches)} match(es)")
    
    for i, match in enumerate(matches[:20], 1):
        with st.expander(f"Match {i}: Position {match.start_1based}"):
            st.markdown(f"**Sequence:** {match.sequence_id}")
            st.markdown(f"**Position:** {match.start_1based} (1-based)")
            st.markdown(f"**Matched:** `{match.matched_sequence}`")
            st.markdown(f"**Context:** `{match.context}`")


def display_stats_badges() -> None:
    """
    Display quick stats badges in the header area.
    """
    st.markdown(f"""
    <div style="display: flex; justify-content: center; gap: 2rem; margin: 1.5rem 0; flex-wrap: wrap;">
        <div style="text-align: center;">
            <div style="color: {COLORS['primary']}; font-size: 1.5rem; font-weight: 600;">10+</div>
            <div style="color: {COLORS['text_secondary']}; font-size: 0.85rem;">Analysis Types</div>
        </div>
        <div style="text-align: center;">
            <div style="color: {COLORS['primary']}; font-size: 1.5rem; font-weight: 600;">2</div>
            <div style="color: {COLORS['text_secondary']}; font-size: 0.85rem;">File Formats</div>
        </div>
        <div style="text-align: center;">
            <div style="color: {COLORS['primary']}; font-size: 1.5rem; font-weight: 600;">⚡</div>
            <div style="color: {COLORS['text_secondary']}; font-size: 0.85rem;">Instant Results</div>
        </div>
    </div>
    """, unsafe_allow_html=True)
