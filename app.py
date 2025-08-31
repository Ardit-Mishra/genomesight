import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO, BytesIO
import zipfile
import json
from datetime import datetime

from sequence_analyzer import SequenceAnalyzer
from visualizations import create_gc_content_plot, create_nucleotide_composition_plot, create_sequence_length_plot, create_quality_scores_plot
from utils import validate_file_format, parse_uploaded_files, generate_report

# Page configuration
st.set_page_config(
    page_title="Genome Sequencing Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.5rem;
        color: #2c3e50;
        margin: 1rem 0;
        border-bottom: 2px solid #1f77b4;
        padding-bottom: 0.5rem;
    }
    .metric-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
</style>
""", unsafe_allow_html=True)

def main():
    # Main header
    st.markdown('<h1 class="main-header">🧬 Genome Sequencing Analyzer</h1>', unsafe_allow_html=True)
    st.markdown("Upload, analyze, and visualize genomic data with comprehensive statistical analysis and downloadable results.")
    
    # Sidebar for file upload and settings
    with st.sidebar:
        st.markdown('<h2 class="section-header">📁 File Upload</h2>', unsafe_allow_html=True)
        
        # File uploader
        uploaded_files = st.file_uploader(
            "Choose genome sequence files",
            type=['fasta', 'fa', 'fastq', 'fq'],
            accept_multiple_files=True,
            help="Upload FASTA (.fasta, .fa) or FASTQ (.fastq, .fq) files"
        )
        
        # Analysis options
        st.markdown('<h3 class="section-header">⚙️ Analysis Options</h3>', unsafe_allow_html=True)
        
        analyze_gc_content = st.checkbox("GC Content Analysis", value=True)
        analyze_composition = st.checkbox("Nucleotide Composition", value=True)
        analyze_quality = st.checkbox("Quality Scores (FASTQ only)", value=True)
        search_pattern = st.text_input("Search Pattern (optional)", placeholder="ATCG")
        
        # Batch processing options
        batch_size = st.selectbox("Batch Processing Size", [10, 50, 100, 500], index=1)
        
    # Main content area
    if uploaded_files:
        # Initialize session state for results
        if 'analysis_results' not in st.session_state:
            st.session_state.analysis_results = {}
        
        # Process files
        with st.spinner("Processing uploaded files..."):
            try:
                processed_files = parse_uploaded_files(uploaded_files)
                
                if processed_files:
                    st.success(f"Successfully loaded {len(processed_files)} files")
                    
                    # Create tabs for different views
                    tab1, tab2, tab3, tab4 = st.tabs(["📊 Overview", "🔬 Detailed Analysis", "📈 Visualizations", "📥 Download Results"])
                    
                    with tab1:
                        display_overview(processed_files)
                    
                    with tab2:
                        display_detailed_analysis(processed_files, analyze_gc_content, analyze_composition, analyze_quality, search_pattern)
                    
                    with tab3:
                        display_visualizations(processed_files, analyze_gc_content, analyze_composition, analyze_quality)
                    
                    with tab4:
                        display_download_section(processed_files)
                
                else:
                    st.error("No valid sequence files found. Please check your file formats.")
                    
            except Exception as e:
                st.error(f"Error processing files: {str(e)}")
                st.info("Please ensure your files are in valid FASTA or FASTQ format.")
    
    else:
        # Welcome screen
        display_welcome_screen()

def display_welcome_screen():
    """Display welcome screen with instructions"""
    st.markdown('<h2 class="section-header">Welcome to Genome Sequencing Analyzer</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        ### 🚀 Features
        - **Multi-format Support**: FASTA and FASTQ files
        - **Comprehensive Analysis**: GC content, nucleotide composition, sequence statistics
        - **Interactive Visualizations**: Charts and plots using Plotly
        - **Batch Processing**: Handle multiple files efficiently
        - **Pattern Search**: Find specific sequences
        - **Downloadable Reports**: Export results in multiple formats
        """)
    
    with col2:
        st.markdown("""
        ### 📋 How to Use
        1. **Upload Files**: Use the sidebar to upload your sequence files
        2. **Select Options**: Choose analysis parameters
        3. **View Results**: Explore overview, detailed analysis, and visualizations
        4. **Download**: Export your results as CSV, JSON, or comprehensive reports
        
        ### 📁 Supported Formats
        - FASTA files (.fasta, .fa)
        - FASTQ files (.fastq, .fq)
        """)
    
    # Sample data section
    st.markdown('<h3 class="section-header">🧪 Example Analysis</h3>', unsafe_allow_html=True)
    st.info("Upload your genome sequence files using the sidebar to begin analysis. The analyzer supports both single and multiple file processing.")

def display_overview(processed_files):
    """Display overview statistics"""
    st.markdown('<h2 class="section-header">📊 File Overview</h2>', unsafe_allow_html=True)
    
    if not processed_files:
        st.warning("No files to analyze")
        return
    
    # Calculate overall statistics
    total_sequences = sum(len(data['sequences']) for data in processed_files.values())
    total_nucleotides = sum(
        sum(len(seq.seq) for seq in data['sequences']) 
        for data in processed_files.values()
    )
    file_types = [data['format'] for data in processed_files.values()]
    
    # Display metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Files Uploaded", len(processed_files))
    
    with col2:
        st.metric("Total Sequences", total_sequences)
    
    with col3:
        st.metric("Total Nucleotides", f"{total_nucleotides:,}")
    
    with col4:
        st.metric("File Types", f"FASTA: {file_types.count('fasta')}, FASTQ: {file_types.count('fastq')}")
    
    # File details table
    st.markdown('<h3 class="section-header">📄 File Details</h3>', unsafe_allow_html=True)
    
    file_details = []
    for filename, data in processed_files.items():
        file_details.append({
            'Filename': filename,
            'Format': data['format'].upper(),
            'Sequences': len(data['sequences']),
            'Total Length': sum(len(seq.seq) for seq in data['sequences']),
            'Avg Length': round(sum(len(seq.seq) for seq in data['sequences']) / len(data['sequences']), 2) if data['sequences'] else 0
        })
    
    df_details = pd.DataFrame(file_details)
    st.dataframe(df_details, use_container_width=True)

def display_detailed_analysis(processed_files, analyze_gc_content, analyze_composition, analyze_quality, search_pattern):
    """Display detailed analysis results"""
    st.markdown('<h2 class="section-header">🔬 Detailed Analysis</h2>', unsafe_allow_html=True)
    
    if not processed_files:
        st.warning("No files to analyze")
        return
    
    # Initialize analyzer
    analyzer = SequenceAnalyzer()
    
    # Process each file
    for filename, data in processed_files.items():
        with st.expander(f"📁 {filename}", expanded=True):
            sequences = data['sequences']
            
            if not sequences:
                st.warning("No sequences found in this file")
                continue
            
            # Basic statistics
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Basic Statistics**")
                total_length = sum(len(seq.seq) for seq in sequences)
                avg_length = total_length / len(sequences) if sequences else 0
                min_length = min(len(seq.seq) for seq in sequences) if sequences else 0
                max_length = max(len(seq.seq) for seq in sequences) if sequences else 0
                
                stats_df = pd.DataFrame({
                    'Metric': ['Sequence Count', 'Total Length', 'Average Length', 'Min Length', 'Max Length'],
                    'Value': [len(sequences), f"{total_length:,}", f"{avg_length:.2f}", min_length, max_length]
                })
                st.dataframe(stats_df, hide_index=True)
            
            with col2:
                if analyze_gc_content:
                    st.markdown("**GC Content Analysis**")
                    gc_contents = [analyzer.calculate_gc_content(seq.seq) for seq in sequences]
                    avg_gc = sum(gc_contents) / len(gc_contents) if gc_contents else 0
                    min_gc = min(gc_contents) if gc_contents else 0
                    max_gc = max(gc_contents) if gc_contents else 0
                    
                    gc_df = pd.DataFrame({
                        'Metric': ['Average GC%', 'Min GC%', 'Max GC%'],
                        'Value': [f"{avg_gc:.2f}%", f"{min_gc:.2f}%", f"{max_gc:.2f}%"]
                    })
                    st.dataframe(gc_df, hide_index=True)
            
            # Nucleotide composition
            if analyze_composition:
                st.markdown("**Nucleotide Composition**")
                combined_sequence = ''.join(str(seq.seq) for seq in sequences)
                composition = analyzer.calculate_nucleotide_composition(combined_sequence)
                
                comp_df = pd.DataFrame(list(composition.items()), columns=['Nucleotide', 'Count'])
                comp_df['Percentage'] = (comp_df['Count'] / comp_df['Count'].sum() * 100).round(2)
                st.dataframe(comp_df, hide_index=True)
            
            # Quality scores for FASTQ
            if analyze_quality and data['format'] == 'fastq':
                st.markdown("**Quality Score Analysis**")
                quality_stats = analyzer.calculate_quality_statistics(sequences)
                if quality_stats:
                    qual_df = pd.DataFrame([quality_stats])
                    st.dataframe(qual_df, hide_index=True)
            
            # Pattern search
            if search_pattern:
                st.markdown(f"**Pattern Search: '{search_pattern}'**")
                matches = analyzer.search_pattern(sequences, search_pattern)
                if matches:
                    match_df = pd.DataFrame(matches)
                    st.dataframe(match_df, hide_index=True)
                else:
                    st.info(f"Pattern '{search_pattern}' not found in any sequences")

def display_visualizations(processed_files, analyze_gc_content, analyze_composition, analyze_quality):
    """Display interactive visualizations"""
    st.markdown('<h2 class="section-header">📈 Visualizations</h2>', unsafe_allow_html=True)
    
    if not processed_files:
        st.warning("No files to visualize")
        return
    
    analyzer = SequenceAnalyzer()
    
    # Collect data for visualizations
    all_sequences = []
    file_labels = []
    
    for filename, data in processed_files.items():
        for seq in data['sequences']:
            all_sequences.append(seq)
            file_labels.append(filename)
    
    if not all_sequences:
        st.warning("No sequences available for visualization")
        return
    
    # Sequence length distribution
    st.markdown("### Sequence Length Distribution")
    lengths = [len(seq.seq) for seq in all_sequences]
    length_fig = create_sequence_length_plot(lengths, file_labels)
    st.plotly_chart(length_fig, use_container_width=True)
    
    # GC content distribution
    if analyze_gc_content:
        st.markdown("### GC Content Distribution")
        gc_contents = [analyzer.calculate_gc_content(seq.seq) for seq in all_sequences]
        gc_fig = create_gc_content_plot(gc_contents, file_labels)
        st.plotly_chart(gc_fig, use_container_width=True)
    
    # Nucleotide composition
    if analyze_composition:
        st.markdown("### Nucleotide Composition by File")
        
        # Create composition data by file
        composition_data = {}
        for filename, data in processed_files.items():
            combined_seq = ''.join(str(seq.seq) for seq in data['sequences'])
            composition_data[filename] = analyzer.calculate_nucleotide_composition(combined_seq)
        
        comp_fig = create_nucleotide_composition_plot(composition_data)
        st.plotly_chart(comp_fig, use_container_width=True)
    
    # Quality scores for FASTQ files
    if analyze_quality:
        fastq_files = {k: v for k, v in processed_files.items() if v['format'] == 'fastq'}
        if fastq_files:
            st.markdown("### Quality Score Distribution (FASTQ files)")
            
            for filename, data in fastq_files.items():
                quality_data = analyzer.get_quality_data(data['sequences'])
                if quality_data:
                    qual_fig = create_quality_scores_plot(quality_data, filename)
                    st.plotly_chart(qual_fig, use_container_width=True)

def display_download_section(processed_files):
    """Display download options"""
    st.markdown('<h2 class="section-header">📥 Download Results</h2>', unsafe_allow_html=True)
    
    if not processed_files:
        st.warning("No results to download")
        return
    
    analyzer = SequenceAnalyzer()
    
    # Generate comprehensive analysis results
    results = {}
    for filename, data in processed_files.items():
        sequences = data['sequences']
        
        if sequences:
            # Basic statistics
            total_length = sum(len(seq.seq) for seq in sequences)
            avg_length = total_length / len(sequences)
            
            # GC content
            gc_contents = [analyzer.calculate_gc_content(seq.seq) for seq in sequences]
            avg_gc = sum(gc_contents) / len(gc_contents)
            
            # Nucleotide composition
            combined_seq = ''.join(str(seq.seq) for seq in sequences)
            composition = analyzer.calculate_nucleotide_composition(combined_seq)
            
            # Quality statistics for FASTQ
            quality_stats = None
            if data['format'] == 'fastq':
                quality_stats = analyzer.calculate_quality_statistics(sequences)
            
            results[filename] = {
                'format': data['format'],
                'sequence_count': len(sequences),
                'total_length': total_length,
                'average_length': avg_length,
                'min_length': min(len(seq.seq) for seq in sequences),
                'max_length': max(len(seq.seq) for seq in sequences),
                'average_gc_content': avg_gc,
                'nucleotide_composition': composition,
                'quality_statistics': quality_stats
            }
    
    # Download options
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("📊 Download CSV Report"):
            csv_data = generate_report(results, 'csv')
            st.download_button(
                label="Download CSV",
                data=csv_data,
                file_name=f"genome_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )
    
    with col2:
        if st.button("📋 Download JSON Report"):
            json_data = generate_report(results, 'json')
            st.download_button(
                label="Download JSON",
                data=json_data,
                file_name=f"genome_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                mime="application/json"
            )
    
    with col3:
        if st.button("📄 Download Full Report"):
            full_report = generate_report(results, 'full')
            st.download_button(
                label="Download Full Report",
                data=full_report,
                file_name=f"genome_analysis_full_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                mime="text/plain"
            )
    
    # Preview results
    st.markdown("### 👀 Results Preview")
    
    # Convert results to DataFrame for display
    preview_data = []
    for filename, file_results in results.items():
        preview_data.append({
            'Filename': filename,
            'Format': file_results['format'].upper(),
            'Sequences': file_results['sequence_count'],
            'Total Length': f"{file_results['total_length']:,}",
            'Avg Length': f"{file_results['average_length']:.2f}",
            'Avg GC%': f"{file_results['average_gc_content']:.2f}%"
        })
    
    if preview_data:
        preview_df = pd.DataFrame(preview_data)
        st.dataframe(preview_df, use_container_width=True)

if __name__ == "__main__":
    main()
