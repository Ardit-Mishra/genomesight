"""
Visualization Module

This module provides functions for creating interactive Plotly charts
for genomic data visualization. All charts use a consistent dark theme
suitable for scientific applications.

Design Notes:
    - All functions return Plotly Figure objects
    - "Laboratory Instrument" dark theme, shared with the rest of the portfolio
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


# Single source of truth for the palette. Charts previously hardcoded their own
# hexes, so the plots and the surrounding UI drifted apart; importing keeps one
# definition. (styles.py has no Streamlit dependency, so this stays importable.)
from app.ui.styles import COLORS, NUCLEOTIDE_COLORS  # noqa: E402

DARK_THEME = {
    'paper_bgcolor': COLORS['background'],
    'plot_bgcolor': COLORS['surface'],
    'font_color': COLORS['text_primary'],
}

#: The accent used for a single data series, and the muted tone for reference
#: lines and thresholds — a threshold is context, not a finding, so it never
#: competes with the data.
SERIES_COLOR = COLORS['primary']
REFERENCE_COLOR = COLORS['text_secondary']
THRESHOLD_COLOR = COLORS['warning']


def apply_dark_theme(fig: go.Figure) -> go.Figure:
    """Apply consistent dark theme to a Plotly figure."""
    fig.update_layout(
        paper_bgcolor=DARK_THEME['paper_bgcolor'],
        plot_bgcolor=DARK_THEME['plot_bgcolor'],
        font=dict(color=DARK_THEME['font_color']),
        xaxis=dict(gridcolor=COLORS['border'], zerolinecolor=COLORS['border']),
        yaxis=dict(gridcolor=COLORS['border'], zerolinecolor=COLORS['border'])
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
        line_color=SERIES_COLOR,
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
        line=dict(color=SERIES_COLOR, width=2),
        fill='tozeroy',
        fillcolor='rgba(110, 155, 255, 0.16)'
    ))
    
    fig.add_hline(
        y=50,
        line_dash="dash",
        line_color=REFERENCE_COLOR,
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

    # Plotly's default autorange leaves just enough headroom for the bar
    # itself, not for a text label sitting above it — on real data the
    # tallest bar's "24.7%" label collided with the plot's top edge and
    # rendered half-clipped. 15% of headroom above the tallest bar gives
    # outside labels somewhere to sit.
    max_pct = max(percentages) if percentages else 0
    fig.update_layout(
        title='Nucleotide Composition',
        xaxis_title='Nucleotide',
        yaxis=dict(range=[0, max_pct * 1.15 if max_pct > 0 else 1]),
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
        marker_color=SERIES_COLOR,
        text=counts,
        textposition='outside'
    ))

    # Same headroom fix as create_composition_plot — an outside text label
    # on the tallest bar otherwise collides with the plot's top edge.
    max_count = max(counts) if counts else 0
    fig.update_layout(
        title=f'Top 20 {k}-mers by Frequency',
        xaxis_title=f'{k}-mer',
        yaxis=dict(range=[0, max_count * 1.15 if max_count > 0 else 1]),
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
        fillcolor='rgba(110, 155, 255, 0.16)',
        line=dict(color='rgba(255,255,255,0)'),
        name='25-75 Percentile'
    ))
    
    fig.add_trace(go.Scatter(
        x=positions,
        y=position_means,
        mode='lines',
        name='Mean Quality',
        line=dict(color=SERIES_COLOR, width=2)
    ))
    
    fig.add_hline(y=20, line_dash="dash", line_color=THRESHOLD_COLOR,
                  annotation_text="Q20")
    fig.add_hline(y=30, line_dash="dash", line_color=THRESHOLD_COLOR,
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

    # nbinsx alone breaks down when every read is the same length (the
    # common case for a fixed-cycle sequencer): with zero span to divide,
    # Plotly still draws 30 bins and lands on a sub-1-bp bin size, so the
    # single true value spreads across fractional bp ticks like "39.6" /
    # "40.2" — precision the integer data never had. An explicit integer
    # bin width, floored at 1 bp, keeps every bin a whole base pair.
    span = max(lengths) - min(lengths) if lengths else 0
    bin_size = max(1, round(span / 30)) if span > 0 else 1
    fig.add_trace(go.Histogram(
        x=lengths,
        xbins=dict(size=bin_size),
        marker_color=SERIES_COLOR,
        opacity=0.7
    ))

    avg_len = np.mean(lengths)
    fig.add_vline(
        x=avg_len,
        line_dash="dash",
        line_color=THRESHOLD_COLOR,
        annotation_text=f"Mean: {avg_len:.0f}"
    )

    fig.update_layout(
        title='Sequence Length Distribution',
        xaxis_title='Sequence Length (bp)',
        yaxis_title='Count',
        height=400
    )

    if span == 0:
        # A single-value dataset still gets autoranged tightly around that
        # value, and Plotly then spaces ticks in fractions of the (tiny)
        # visible span — "39.6", "40.2" bp. Pad the range and pin the tick
        # step to 1 bp so the axis reads as the whole numbers it actually is.
        fig.update_xaxes(range=[lengths[0] - 3, lengths[0] + 3], dtick=1)

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
    
    # Same fixed-bin-width fix as create_length_distribution: nbinsx alone
    # produces fractional-bp bins when every ORF happens to be the same
    # length, which is exactly the "min_length" default filtering can
    # produce on a short, low-diversity input.
    lengths = [orf.length_nt for orf in orfs]
    span = max(lengths) - min(lengths) if lengths else 0
    bin_size = max(1, round(span / 20)) if span > 0 else 1
    fig.add_trace(
        go.Histogram(x=lengths, xbins=dict(size=bin_size), marker_color=SERIES_COLOR),
        row=1, col=1
    )
    
    frames = Counter([f"{orf.strand}{orf.frame}" for orf in orfs])
    # Fixed, biologically ordered set of the six reading frames rather than
    # whatever order Counter first saw them in — and only the frames that
    # actually occur, so a forward-only run doesn't show three empty '-' bars.
    frame_order = [f"{s}{f}" for s in ('+', '-') for f in (1, 2, 3)]
    frame_labels = [f for f in frame_order if f in frames]
    frame_values = [frames[f] for f in frame_labels]
    fig.add_trace(
        go.Bar(x=frame_labels, y=frame_values,
               marker_color=SERIES_COLOR),
        row=1, col=2
    )
    # Labels like "+1" / "-1" parse as numbers, so Plotly's default axis-type
    # inference plots them on a continuous number line — "+1" and "-1" land
    # symmetrically around a "0" that doesn't correspond to any frame, and
    # adjacent frames end up spaced apart by magnitude instead of shown as
    # the six discrete bins they are. Force category so each frame is its
    # own even bar.
    fig.update_xaxes(type='category', row=1, col=2)
    if span == 0:
        # Same tight-autorange fractional-tick artifact as
        # create_length_distribution, for the same reason: one bin, zero
        # span, Plotly zooms in and ticks it in fractions of a base pair.
        fig.update_xaxes(range=[lengths[0] - 3, lengths[0] + 3], dtick=1, row=1, col=1)

    fig.update_layout(
        title='ORF Analysis',
        height=400,
        showlegend=False
    )

    return apply_dark_theme(fig)
