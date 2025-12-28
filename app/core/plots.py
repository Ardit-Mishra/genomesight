"""
Visualization Module

This module provides functions for creating interactive Plotly charts
for genomic data visualization. All charts use a consistent dark theme
suitable for scientific applications.

Design Notes:
    - All functions return Plotly Figure objects
    - Dark theme with teal accent colors
    - Interactive by default with hover information
    - No Streamlit dependencies

Chart Types:
    - GC content distribution and sliding window
    - Nucleotide composition bar charts
    - K-mer frequency heatmaps
    - Quality score distributions
    - Sequence length distributions
"""

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from typing import List, Dict, Any, Optional
from collections import Counter
import logging

logger = logging.getLogger(__name__)


DARK_THEME = {
    'paper_bgcolor': 'rgba(26, 29, 41, 1)',
    'plot_bgcolor': 'rgba(37, 41, 60, 1)',
    'font_color': '#E8EAED'
}

NUCLEOTIDE_COLORS = {
    'A': '#FF6B9D',
    'T': '#FFC252',
    'C': '#00C9A7',
    'G': '#4D96FF',
    'N': '#9CA3AF'
}


def apply_dark_theme(fig: go.Figure) -> go.Figure:
    """Apply consistent dark theme to a Plotly figure."""
    fig.update_layout(
        paper_bgcolor=DARK_THEME['paper_bgcolor'],
        plot_bgcolor=DARK_THEME['plot_bgcolor'],
        font=dict(color=DARK_THEME['font_color']),
        xaxis=dict(gridcolor='#373B4D', zerolinecolor='#373B4D'),
        yaxis=dict(gridcolor='#373B4D', zerolinecolor='#373B4D')
    )
    return fig


def create_gc_plot(gc_contents: List[float], labels: List[str]) -> go.Figure:
    """
    Create GC content distribution histogram.
    
    Args:
        gc_contents: List of GC percentages
        labels: List of labels for each value
    
    Returns:
        Plotly Figure object
    """
    df = pd.DataFrame({'GC_Content': gc_contents, 'Label': labels})
    
    fig = px.histogram(
        df,
        x='GC_Content',
        color='Label',
        title='GC Content Distribution',
        labels={'GC_Content': 'GC Content (%)', 'count': 'Count'},
        nbins=30,
        opacity=0.7
    )
    
    avg_gc = np.mean(gc_contents)
    fig.add_vline(
        x=avg_gc,
        line_dash="dash",
        line_color="#00C9A7",
        annotation_text=f"Average: {avg_gc:.1f}%"
    )
    
    fig.update_layout(
        height=400,
        xaxis_title="GC Content (%)",
        yaxis_title="Number of Sequences"
    )
    
    return apply_dark_theme(fig)


def create_gc_sliding_window(
    sequence: str,
    window_size: int = 100,
    step: int = 10
) -> go.Figure:
    """
    Create sliding window GC content plot.
    
    Args:
        sequence: DNA sequence string
        window_size: Window size in bases
        step: Step size between windows
    
    Returns:
        Plotly Figure object
    """
    positions = []
    gc_values = []
    
    for i in range(0, len(sequence) - window_size + 1, step):
        window = sequence[i:i + window_size].upper()
        gc = (window.count('G') + window.count('C'))
        total = len([b for b in window if b in 'ATCG'])
        gc_pct = (gc / total * 100) if total > 0 else 0
        
        positions.append(i + window_size // 2)
        gc_values.append(gc_pct)
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=positions,
        y=gc_values,
        mode='lines',
        name='GC Content',
        line=dict(color='#00C9A7', width=2),
        fill='tozeroy',
        fillcolor='rgba(0, 201, 167, 0.2)'
    ))
    
    fig.add_hline(
        y=50,
        line_dash="dash",
        line_color="#9CA3AF",
        annotation_text="50%"
    )
    
    fig.update_layout(
        title=f'GC Content (Window: {window_size} bp)',
        xaxis_title='Position (bp)',
        yaxis_title='GC Content (%)',
        height=400
    )
    
    return apply_dark_theme(fig)


def create_composition_plot(composition: Dict[str, int]) -> go.Figure:
    """
    Create nucleotide composition bar chart.
    
    Args:
        composition: Dictionary of nucleotide counts
    
    Returns:
        Plotly Figure object
    """
    nucleotides = ['A', 'T', 'C', 'G', 'N']
    counts = [composition.get(n, 0) for n in nucleotides]
    total = sum(counts)
    percentages = [(c / total * 100) if total > 0 else 0 for c in counts]
    colors = [NUCLEOTIDE_COLORS[n] for n in nucleotides]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=nucleotides,
        y=percentages,
        marker_color=colors,
        text=[f'{p:.1f}%' for p in percentages],
        textposition='outside'
    ))
    
    fig.update_layout(
        title='Nucleotide Composition',
        xaxis_title='Nucleotide',
        yaxis_title='Percentage (%)',
        height=400,
        showlegend=False
    )
    
    return apply_dark_theme(fig)


def create_kmer_plot(kmer_counts: Dict[str, int], k: int) -> go.Figure:
    """
    Create k-mer frequency bar chart.
    
    Args:
        kmer_counts: Dictionary of k-mer counts
        k: K-mer size
    
    Returns:
        Plotly Figure object
    """
    sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)[:20]
    kmers = [x[0] for x in sorted_kmers]
    counts = [x[1] for x in sorted_kmers]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=kmers,
        y=counts,
        marker_color='#00C9A7',
        text=counts,
        textposition='outside'
    ))
    
    fig.update_layout(
        title=f'Top 20 {k}-mers by Frequency',
        xaxis_title=f'{k}-mer',
        yaxis_title='Count',
        height=400,
        xaxis_tickangle=-45
    )
    
    return apply_dark_theme(fig)


def create_quality_plot(quality_data: List[List[int]]) -> go.Figure:
    """
    Create quality score distribution plot for FASTQ data.
    
    Args:
        quality_data: List of quality scores per sequence
    
    Returns:
        Plotly Figure object
    """
    if not quality_data:
        return go.Figure()
    
    max_length = max(len(q) for q in quality_data)
    
    position_means = []
    position_q25 = []
    position_q75 = []
    positions = []
    
    for pos in range(max_length):
        scores = [q[pos] for q in quality_data if pos < len(q)]
        if scores:
            positions.append(pos + 1)
            position_means.append(np.mean(scores))
            position_q25.append(np.percentile(scores, 25))
            position_q75.append(np.percentile(scores, 75))
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=positions + positions[::-1],
        y=position_q75 + position_q25[::-1],
        fill='toself',
        fillcolor='rgba(0, 201, 167, 0.2)',
        line=dict(color='rgba(255,255,255,0)'),
        name='25-75 Percentile'
    ))
    
    fig.add_trace(go.Scatter(
        x=positions,
        y=position_means,
        mode='lines',
        name='Mean Quality',
        line=dict(color='#00C9A7', width=2)
    ))
    
    fig.add_hline(y=20, line_dash="dash", line_color="#FFC252",
                  annotation_text="Q20")
    fig.add_hline(y=30, line_dash="dash", line_color="#FF6B9D",
                  annotation_text="Q30")
    
    fig.update_layout(
        title='Quality Scores by Position',
        xaxis_title='Position',
        yaxis_title='Quality Score',
        height=400
    )
    
    return apply_dark_theme(fig)


def create_length_distribution(lengths: List[int]) -> go.Figure:
    """
    Create sequence length distribution histogram.
    
    Args:
        lengths: List of sequence lengths
    
    Returns:
        Plotly Figure object
    """
    fig = go.Figure()
    
    fig.add_trace(go.Histogram(
        x=lengths,
        nbinsx=30,
        marker_color='#00C9A7',
        opacity=0.7
    ))
    
    avg_len = np.mean(lengths)
    fig.add_vline(
        x=avg_len,
        line_dash="dash",
        line_color="#FFC252",
        annotation_text=f"Mean: {avg_len:.0f}"
    )
    
    fig.update_layout(
        title='Sequence Length Distribution',
        xaxis_title='Sequence Length (bp)',
        yaxis_title='Count',
        height=400
    )
    
    return apply_dark_theme(fig)


def create_orf_plot(orfs: List) -> go.Figure:
    """
    Create ORF visualization plot.
    
    Args:
        orfs: List of ORF objects
    
    Returns:
        Plotly Figure object
    """
    if not orfs:
        return go.Figure()
    
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=('ORF Length Distribution', 'ORFs by Frame'))
    
    lengths = [orf.length_nt for orf in orfs]
    fig.add_trace(
        go.Histogram(x=lengths, nbinsx=20, marker_color='#00C9A7'),
        row=1, col=1
    )
    
    frames = Counter([f"{orf.strand}{orf.frame}" for orf in orfs])
    fig.add_trace(
        go.Bar(x=list(frames.keys()), y=list(frames.values()),
               marker_color='#4D96FF'),
        row=1, col=2
    )
    
    fig.update_layout(
        title='ORF Analysis',
        height=400,
        showlegend=False
    )
    
    return apply_dark_theme(fig)
