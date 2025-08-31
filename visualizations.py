import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
from typing import List, Dict, Any
import streamlit as st

def create_gc_content_plot(gc_contents: List[float], file_labels: List[str]) -> go.Figure:
    """
    Create GC content distribution plot
    
    Args:
        gc_contents (List[float]): List of GC content percentages
        file_labels (List[str]): List of file labels for each sequence
        
    Returns:
        go.Figure: Plotly figure
    """
    df = pd.DataFrame({
        'GC_Content': gc_contents,
        'File': file_labels
    })
    
    fig = px.histogram(
        df, 
        x='GC_Content', 
        color='File',
        title='GC Content Distribution',
        labels={'GC_Content': 'GC Content (%)', 'count': 'Number of Sequences'},
        nbins=30,
        opacity=0.7
    )
    
    # Add average line
    avg_gc = np.mean(gc_contents)
    fig.add_vline(
        x=avg_gc, 
        line_dash="dash", 
        line_color="red",
        annotation_text=f"Average: {avg_gc:.2f}%"
    )
    
    fig.update_layout(
        xaxis_title="GC Content (%)",
        yaxis_title="Number of Sequences",
        showlegend=True,
        height=500,
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#fafafa")
    )
    
    return fig

def create_nucleotide_composition_plot(composition_data: Dict[str, Dict[str, int]]) -> go.Figure:
    """
    Create nucleotide composition plot
    
    Args:
        composition_data (Dict[str, Dict[str, int]]): Composition data by file
        
    Returns:
        go.Figure: Plotly figure
    """
    # Prepare data for stacked bar chart
    files = list(composition_data.keys())
    nucleotides = ['A', 'T', 'C', 'G', 'N']
    
    fig = go.Figure()
    
    colors = {
        'A': '#FF6B6B',  # Red
        'T': '#4ECDC4',  # Teal
        'C': '#45B7D1',  # Blue
        'G': '#96CEB4',  # Green
        'N': '#FECA57'   # Yellow
    }
    
    for nucleotide in nucleotides:
        values = []
        for file in files:
            total = sum(composition_data[file].values())
            percentage = (composition_data[file].get(nucleotide, 0) / total * 100) if total > 0 else 0
            values.append(percentage)
        
        fig.add_trace(go.Bar(
            name=nucleotide,
            x=files,
            y=values,
            marker_color=colors.get(nucleotide, '#95A5A6')
        ))
    
    fig.update_layout(
        title='Nucleotide Composition by File',
        xaxis_title='File',
        yaxis_title='Percentage (%)',
        barmode='stack',
        height=500,
        showlegend=True,
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#fafafa")
    )
    
    return fig

def create_sequence_length_plot(lengths: List[int], file_labels: List[str]) -> go.Figure:
    """
    Create sequence length distribution plot
    
    Args:
        lengths (List[int]): List of sequence lengths
        file_labels (List[str]): List of file labels for each sequence
        
    Returns:
        go.Figure: Plotly figure
    """
    df = pd.DataFrame({
        'Length': lengths,
        'File': file_labels
    })
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('Sequence Length Distribution', 'Box Plot by File'),
        vertical_spacing=0.12,
        specs=[[{"secondary_y": False}], [{"secondary_y": False}]]
    )
    
    # Histogram
    for file in df['File'].unique():
        file_data = df[df['File'] == file]
        fig.add_trace(
            go.Histogram(
                x=file_data['Length'],
                name=file,
                opacity=0.7,
                nbinsx=30
            ),
            row=1, col=1
        )
    
    # Box plot
    fig.add_trace(
        go.Box(
            x=df['File'],
            y=df['Length'],
            name='Length Distribution',
            showlegend=False
        ),
        row=2, col=1
    )
    
    fig.update_layout(
        title='Sequence Length Analysis',
        height=800,
        showlegend=True,
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#fafafa")
    )
    
    fig.update_xaxes(title_text="Sequence Length (bp)", row=1, col=1)
    fig.update_yaxes(title_text="Number of Sequences", row=1, col=1)
    fig.update_xaxes(title_text="File", row=2, col=1)
    fig.update_yaxes(title_text="Sequence Length (bp)", row=2, col=1)
    
    return fig

def create_quality_scores_plot(quality_data: List[List[float]], filename: str) -> go.Figure:
    """
    Create quality scores plot for FASTQ data
    
    Args:
        quality_data (List[List[float]]): Quality scores per sequence
        filename (str): Filename for the plot title
        
    Returns:
        go.Figure: Plotly figure
    """
    if not quality_data:
        return go.Figure()
    
    # Calculate position-wise statistics
    max_length = max(len(seq_qual) for seq_qual in quality_data)
    position_stats = []
    
    for pos in range(max_length):
        scores_at_pos = []
        for seq_qual in quality_data:
            if pos < len(seq_qual):
                scores_at_pos.append(seq_qual[pos])
        
        if scores_at_pos:
            position_stats.append({
                'position': pos + 1,
                'mean': np.mean(scores_at_pos),
                'q25': np.percentile(scores_at_pos, 25),
                'median': np.median(scores_at_pos),
                'q75': np.percentile(scores_at_pos, 75),
                'min': min(scores_at_pos),
                'max': max(scores_at_pos)
            })
    
    df_stats = pd.DataFrame(position_stats)
    
    fig = go.Figure()
    
    # Add mean line
    fig.add_trace(go.Scatter(
        x=df_stats['position'],
        y=df_stats['mean'],
        mode='lines',
        name='Mean Quality',
        line=dict(color='blue', width=2)
    ))
    
    # Add median line
    fig.add_trace(go.Scatter(
        x=df_stats['position'],
        y=df_stats['median'],
        mode='lines',
        name='Median Quality',
        line=dict(color='green', width=2)
    ))
    
    # Add quartile range
    fig.add_trace(go.Scatter(
        x=df_stats['position'].tolist() + df_stats['position'].tolist()[::-1],
        y=df_stats['q75'].tolist() + df_stats['q25'].tolist()[::-1],
        fill='toself',
        fillcolor='rgba(0,100,80,0.2)',
        line=dict(color='rgba(255,255,255,0)'),
        name='25-75 Percentile',
        showlegend=True
    ))
    
    # Add quality threshold lines
    fig.add_hline(y=20, line_dash="dash", line_color="orange", 
                  annotation_text="Q20 threshold")
    fig.add_hline(y=30, line_dash="dash", line_color="red", 
                  annotation_text="Q30 threshold")
    
    fig.update_layout(
        title=f'Quality Scores by Position - {filename}',
        xaxis_title='Position in Read',
        yaxis_title='Quality Score',
        height=500,
        showlegend=True,
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#fafafa")
    )
    
    return fig

def create_gc_content_heatmap(sequences: List, window_size: int = 100) -> go.Figure:
    """
    Create GC content heatmap along sequences
    
    Args:
        sequences (List): List of Bio.SeqRecord objects
        window_size (int): Window size for GC calculation
        
    Returns:
        go.Figure: Plotly figure
    """
    from sequence_analyzer import SequenceAnalyzer
    analyzer = SequenceAnalyzer()
    
    gc_matrix = []
    sequence_names = []
    
    for i, seq in enumerate(sequences[:50]):  # Limit to first 50 sequences
        sequence_str = str(seq.seq)
        seq_name = seq.id if hasattr(seq, 'id') else f"Sequence_{i}"
        sequence_names.append(seq_name)
        
        gc_values = []
        for j in range(0, len(sequence_str) - window_size + 1, window_size):
            window = sequence_str[j:j + window_size]
            gc_content = analyzer.calculate_gc_content(window)
            gc_values.append(gc_content)
        
        gc_matrix.append(gc_values)
    
    if not gc_matrix:
        return go.Figure()
    
    # Pad shorter sequences with NaN
    max_windows = max(len(row) for row in gc_matrix)
    for row in gc_matrix:
        while len(row) < max_windows:
            row.append(np.nan)
    
    fig = go.Figure(data=go.Heatmap(
        z=gc_matrix,
        x=[f"Window {i+1}" for i in range(max_windows)],
        y=sequence_names,
        colorscale='RdYlBu_r',
        hoverongaps=False,
        colorbar=dict(title="GC Content (%)")
    ))
    
    fig.update_layout(
        title=f'GC Content Heatmap (Window size: {window_size} bp)',
        xaxis_title='Genomic Windows',
        yaxis_title='Sequences',
        height=600
    )
    
    return fig

def create_pattern_distribution_plot(pattern_matches: List[Dict[str, Any]]) -> go.Figure:
    """
    Create pattern match distribution plot
    
    Args:
        pattern_matches (List[Dict[str, Any]]): Pattern match results
        
    Returns:
        go.Figure: Plotly figure
    """
    if not pattern_matches:
        return go.Figure()
    
    df_matches = pd.DataFrame(pattern_matches)
    
    # Count matches per sequence
    match_counts = df_matches.groupby('sequence_id').size().reset_index(name='match_count')
    
    fig = px.bar(
        match_counts,
        x='sequence_id',
        y='match_count',
        title='Pattern Matches per Sequence',
        labels={'sequence_id': 'Sequence ID', 'match_count': 'Number of Matches'}
    )
    
    fig.update_layout(
        xaxis_title="Sequence ID",
        yaxis_title="Number of Matches",
        height=500,
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#fafafa")
    )
    
    return fig

def create_interactive_sequence_logo(sequences: List, title: str = "Sequence Logo") -> go.Figure:
    """
    Create an interactive sequence logo visualization
    
    Args:
        sequences: List of Bio.SeqRecord objects
        title: Plot title
        
    Returns:
        go.Figure: Plotly figure
    """
    if not sequences:
        return go.Figure()
    
    # Convert sequences to strings and find common length
    seq_strings = [str(seq.seq) for seq in sequences[:100]]  # Limit for performance
    min_length = min(len(seq) for seq in seq_strings)
    
    if min_length == 0:
        return go.Figure()
    
    # Calculate position-wise nucleotide frequencies
    position_data = []
    nucleotides = ['A', 'T', 'C', 'G']
    
    for pos in range(min_length):
        pos_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        valid_count = 0
        
        for seq in seq_strings:
            if pos < len(seq) and seq[pos].upper() in nucleotides:
                pos_counts[seq[pos].upper()] += 1
                valid_count += 1
        
        if valid_count > 0:
            for nuc in nucleotides:
                freq = pos_counts[nuc] / valid_count
                position_data.append({
                    'Position': pos + 1,
                    'Nucleotide': nuc,
                    'Frequency': freq,
                    'Count': pos_counts[nuc]
                })
    
    df = pd.DataFrame(position_data)
    
    # Create stacked bar chart as sequence logo approximation
    fig = px.bar(
        df, 
        x='Position', 
        y='Frequency', 
        color='Nucleotide',
        title=title,
        color_discrete_map={
            'A': '#FF6B6B',  # Red
            'T': '#4ECDC4',  # Teal  
            'C': '#45B7D1',  # Blue
            'G': '#96CEB4'   # Green
        }
    )
    
    fig.update_layout(
        xaxis_title="Position",
        yaxis_title="Nucleotide Frequency",
        barmode='stack',
        height=400,
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#fafafa"),
        showlegend=True
    )
    
    return fig

def create_comparative_analysis_plot(file_data: Dict[str, Dict]) -> go.Figure:
    """
    Create comparative analysis across multiple files
    
    Args:
        file_data: Dictionary of file analysis results
        
    Returns:
        go.Figure: Plotly figure
    """
    if not file_data:
        return go.Figure()
    
    # Prepare comparison data
    comparison_data = []
    
    for filename, data in file_data.items():
        if 'sequences' in data:
            sequences = data['sequences']
            lengths = [len(seq.seq) for seq in sequences]
            
            comparison_data.append({
                'File': filename,
                'Sequence_Count': len(sequences),
                'Avg_Length': np.mean(lengths) if lengths else 0,
                'Total_Length': sum(lengths),
                'Min_Length': min(lengths) if lengths else 0,
                'Max_Length': max(lengths) if lengths else 0
            })
    
    df = pd.DataFrame(comparison_data)
    
    # Create subplots for comparison
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Sequence Count', 'Average Length', 'Total Length', 'Length Range'),
        specs=[[{"type": "bar"}, {"type": "bar"}],
               [{"type": "bar"}, {"type": "scatter"}]]
    )
    
    # Sequence count
    fig.add_trace(
        go.Bar(x=df['File'], y=df['Sequence_Count'], name='Count', marker_color='#00d4aa'),
        row=1, col=1
    )
    
    # Average length
    fig.add_trace(
        go.Bar(x=df['File'], y=df['Avg_Length'], name='Avg Length', marker_color='#45B7D1'),
        row=1, col=2
    )
    
    # Total length
    fig.add_trace(
        go.Bar(x=df['File'], y=df['Total_Length'], name='Total Length', marker_color='#96CEB4'),
        row=2, col=1
    )
    
    # Length range (min-max)
    fig.add_trace(
        go.Scatter(
            x=df['File'], 
            y=df['Min_Length'], 
            mode='markers', 
            name='Min Length',
            marker=dict(color='#FF6B6B', size=10)
        ),
        row=2, col=2
    )
    
    fig.add_trace(
        go.Scatter(
            x=df['File'], 
            y=df['Max_Length'], 
            mode='markers', 
            name='Max Length',
            marker=dict(color='#FECA57', size=10)
        ),
        row=2, col=2
    )
    
    fig.update_layout(
        height=800,
        title='File Comparison Dashboard',
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#fafafa"),
        showlegend=False
    )
    
    return fig
