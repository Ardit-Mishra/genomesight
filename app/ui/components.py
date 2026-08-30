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
from .styles import COLORS, NUCLEOTIDE_COLORS, FONT_MONO, format_sequence_html
from .icons import icon, inline


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
    # Was a 2.5rem gradient panel with a 2px dashed accent border, a teal glow
    # and a 3rem DNA EMOJI over a 1.5rem accent heading — a drop zone shouting
    # louder than the results it exists to produce, sitting directly above
    # Streamlit's own uploader (so the page showed two drop zones in a row).
    # Now a single label line; the real uploader beneath it does the work.
    st.markdown(
        f"""
        <div style="display: flex; align-items: center; gap: 0.6rem; margin: 0 0 0.5rem 0;">
            {icon('sequenceFile', 15, COLORS['text_secondary'])}
            <span style="font-family: {FONT_MONO}; font-size: 0.68rem; letter-spacing: 0.12em;
                         text-transform: uppercase; color: {COLORS['text_secondary']};">Sequence file</span>
            <span style="font-size: 0.8rem; color: {COLORS['text_secondary']};">
                FASTA (.fasta, .fa, .fna) &middot; FASTQ (.fastq, .fq) &middot; GenBank (.gb, .gbk, .genbank)
            </span>
        </div>
        """,
        unsafe_allow_html=True,
    )


def render_empty_state() -> None:
    """Teach the interface before any file is loaded.

    Shown only pre-upload — ``main.py`` gates this on
    ``st.session_state.sequences is None``. Once real results exist below,
    repeating "here's what this tool does" under them would just be clutter,
    so this panel disappears the moment a file (or the sample set) loads.
    """
    formats = "".join(
        f'<div style="margin:.4rem 0">{inline("sequenceFile", text)}</div>'
        for text in (
            "FASTA — .fasta .fa .fna",
            "FASTQ — .fastq .fq, includes per-base quality scores",
            "GenBank — .gb .gbk .genbank, includes feature annotations",
        )
    )
    analyses = "".join(
        f'<div style="margin:.4rem 0">{inline(name, text)}</div>'
        for name, text in (
            ("composition", "GC content & base composition"),
            ("kmer", "k-mer frequency — sidebar, adjustable size"),
            ("readingFrame", "Open reading frame detection — sidebar"),
            ("readingFrame", "Codon usage tables & RSCU — sidebar"),
            ("motif", "Motif & restriction-site search — sidebar"),
            ("alignment", "Pairwise alignment — needs 2+ sequences"),
        )
    )
    st.markdown(
        f"""
        <div style="border: 1px solid {COLORS['border']}; border-radius: 6px;
                    padding: 1.25rem 1.4rem; margin: 0 0 1.25rem 0;
                    display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
                    gap: 1.5rem 2.5rem;">
            <div>
                <div style="display: flex; align-items: center; gap: 0.55rem; margin-bottom: 0.6rem;">
                    {icon('sequenceFile', 20, COLORS['text_primary'])}
                    <span style="font-weight: 600; color: {COLORS['text_primary']};">Accepted formats</span>
                </div>
                {formats}
            </div>
            <div>
                <div style="display: flex; align-items: center; gap: 0.55rem; margin-bottom: 0.6rem;">
                    {icon('ruler', 20, COLORS['text_primary'])}
                    <span style="font-weight: 600; color: {COLORS['text_primary']};">Runs on upload</span>
                </div>
                {analyses}
                <div style="margin-top: 0.7rem; font-size: 0.82rem; color: {COLORS['text_secondary']};">
                    No file on hand? Use "Try Sample File" below to explore with
                    real example sequences instead.
                </div>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def sidebar_section_header(name: str, title: str) -> None:
    """A sidebar tool-group heading with its domain glyph.

    Renders as a literal ``<h4>`` so it inherits the sidebar's mono,
    uppercase, letter-spaced "instrument label" styling from ``styles.py`` —
    a plain styled ``<span>`` would not pick up the
    ``section[data-testid="stSidebar"] h4`` rule.
    """
    st.markdown(
        f'<h4 style="display:flex;align-items:center;gap:.5rem;margin:1rem 0 .6rem 0">'
        f'{icon(name, 15)}<span>{title}</span></h4>',
        unsafe_allow_html=True,
    )


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
    Capability strip under the masthead.

    Was three centred badges reading "10+ / Analysis Types", "2 / File Formats"
    and a lightning EMOJI over "Instant Results". Problems, in order: "10+" is a
    number nobody counted, "Instant Results" is marketing rather than a fact
    about the tool, an emoji is not an icon, and centring them fought the
    left-aligned masthead directly above.

    Now: left-aligned, monospace, and every item is checkable — the analyses are
    named, the formats are named, and "deterministic" is a property of the code
    rather than a claim about speed.
    """
    items = [
        ("Analyses", "GC · composition · k-mer · ORF · codon usage · motif · alignment"),
        ("Formats", "FASTA · FASTQ · GenBank"),
        ("Output", "deterministic — same input, same result"),
    ]
    cells = "".join(
        f"""
        <div style="min-width: 0;">
            <div style="font-family: {FONT_MONO}; font-size: 0.68rem; letter-spacing: 0.12em;
                        text-transform: uppercase; color: {COLORS['text_secondary']};">{label}</div>
            <div style="font-size: 0.9rem; color: {COLORS['text_primary']}; margin-top: 0.25rem;">{value}</div>
        </div>"""
        for label, value in items
    )
    st.markdown(
        f"""
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                    gap: 1.5rem 2.5rem; margin: 0 0 2rem 0; padding-bottom: 1.5rem;
                    border-bottom: 1px solid {COLORS['border']};">{cells}
        </div>
        """,
        unsafe_allow_html=True,
    )
