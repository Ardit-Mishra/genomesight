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

# Custom CSS for better dark theme styling
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #00d4aa;
        text-align: center;
        margin-bottom: 2rem;
        text-shadow: 0 0 10px rgba(0, 212, 170, 0.3);
        animation: glow 2s ease-in-out infinite alternate;
    }
    
    @keyframes glow {
        from { text-shadow: 0 0 10px rgba(0, 212, 170, 0.3); }
        to { text-shadow: 0 0 20px rgba(0, 212, 170, 0.6); }
    }
    
    .section-header {
        font-size: 1.5rem;
        color: #00d4aa;
        margin: 1rem 0;
        border-bottom: 2px solid #00d4aa;
        padding-bottom: 0.5rem;
        background: linear-gradient(90deg, #00d4aa, transparent);
        background-size: 100% 3px;
        background-repeat: no-repeat;
        background-position: bottom;
    }
    
    .metric-card {
        background: linear-gradient(135deg, #262730, #1e2130);
        padding: 1.5rem;
        border-radius: 1rem;
        border: 1px solid #00d4aa;
        box-shadow: 0 4px 15px rgba(0, 212, 170, 0.1);
        transition: transform 0.3s ease, box-shadow 0.3s ease;
    }
    
    .metric-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 8px 25px rgba(0, 212, 170, 0.2);
    }
    
    .stFileUploader > div > div > div {
        background: linear-gradient(135deg, #262730, #1e2130);
        border: 2px dashed #00d4aa;
        border-radius: 1rem;
        padding: 2rem;
        transition: all 0.3s ease;
    }
    
    .stFileUploader > div > div > div:hover {
        border-color: #00ff94;
        background: linear-gradient(135deg, #2a2f3f, #232741);
    }
    
    .stSelectbox > div > div {
        background-color: #262730;
        border: 1px solid #00d4aa;
        border-radius: 0.5rem;
    }
    
    .stCheckbox > label {
        color: #fafafa;
        font-weight: 500;
    }
    
    .stTextInput > div > div > input {
        background-color: #262730;
        border: 1px solid #00d4aa;
        border-radius: 0.5rem;
        color: #fafafa;
    }
    
    .stButton > button {
        background: linear-gradient(135deg, #00d4aa, #00b894);
        color: white;
        border: none;
        border-radius: 0.5rem;
        padding: 0.5rem 1rem;
        font-weight: 600;
        transition: all 0.3s ease;
        box-shadow: 0 4px 15px rgba(0, 212, 170, 0.3);
    }
    
    .stButton > button:hover {
        background: linear-gradient(135deg, #00ff94, #00d4aa);
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(0, 212, 170, 0.4);
    }
    
    .stExpander {
        background: linear-gradient(135deg, #262730, #1e2130);
        border: 1px solid #00d4aa;
        border-radius: 1rem;
        margin: 1rem 0;
    }
    
    .uploadedFile {
        background: linear-gradient(135deg, #262730, #1e2130);
        border: 1px solid #00d4aa;
        border-radius: 0.5rem;
        padding: 0.5rem;
        margin: 0.25rem 0;
    }
    
    .stDataFrame {
        background-color: #262730;
        border-radius: 0.5rem;
        overflow: hidden;
    }
    
    .stAlert {
        background: linear-gradient(135deg, #262730, #1e2130);
        border-left: 4px solid #00d4aa;
        border-radius: 0.5rem;
    }
    
    /* Progress bars */
    .stProgress > div > div > div {
        background: linear-gradient(90deg, #00d4aa, #00ff94);
    }
    
    /* Sidebar styling */
    .css-1d391kg {
        background: linear-gradient(180deg, #1e2130, #262730);
    }
    
    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
    
    .stTabs [data-baseweb="tab"] {
        background: linear-gradient(135deg, #262730, #1e2130);
        border: 1px solid #00d4aa;
        border-radius: 0.5rem;
        padding: 0.5rem 1rem;
        color: #fafafa;
    }
    
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #00d4aa, #00b894);
        color: white;
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
        
        # Interactive file upload guidance
        st.info("💡 **Tip**: Drag and drop your files or click to browse")
        
        # File uploader with enhanced help
        uploaded_files = st.file_uploader(
            "Choose genome sequence files",
            type=['fasta', 'fa', 'fastq', 'fq'],
            accept_multiple_files=True,
            help="📋 **Supported formats:**\n• FASTA (.fasta, .fa)\n• FASTQ (.fastq, .fq)\n\n🔢 **Multiple files:** Upload up to 50 files at once for batch analysis"
        )
        
        # File upload status
        if uploaded_files:
            st.success(f"✅ {len(uploaded_files)} file(s) uploaded successfully!")
            with st.expander("📋 Uploaded Files", expanded=False):
                for file in uploaded_files:
                    file_size = len(file.getvalue()) / 1024  # KB
                    st.write(f"📄 **{file.name}** ({file_size:.1f} KB)")
        else:
            st.warning("⏳ No files uploaded yet")
        
        # Analysis options
        st.markdown('<h3 class="section-header">⚙️ Analysis Options</h3>', unsafe_allow_html=True)
        
        st.markdown("**🔬 Core Analysis**")
        analyze_gc_content = st.checkbox(
            "GC Content Analysis", 
            value=True,
            help="Calculate GC content percentage for each sequence - important for genome characterization"
        )
        analyze_composition = st.checkbox(
            "Nucleotide Composition", 
            value=True,
            help="Count and visualize A, T, C, G, N nucleotides in your sequences"
        )
        analyze_quality = st.checkbox(
            "Quality Scores (FASTQ only)", 
            value=True,
            help="Analyze Phred quality scores for FASTQ files to assess sequencing quality"
        )
        
        st.markdown("**🔍 Pattern Search**")
        search_pattern = st.text_input(
            "Search DNA Pattern", 
            placeholder="e.g., ATCG, TATA, or any sequence",
            help="Search for specific DNA sequences or motifs across all uploaded files"
        )
        
        if search_pattern:
            st.info(f"🔍 Will search for pattern: **{search_pattern.upper()}**")
        
        # Advanced options
        st.markdown("**⚡ Processing Options**")
        batch_size = st.selectbox(
            "Batch Processing Size", 
            [10, 50, 100, 500], 
            index=1,
            help="Number of sequences to process at once - larger batches are faster but use more memory"
        )
        
        # Analysis button
        if uploaded_files:
            if st.button("🚀 Start Analysis", type="primary", use_container_width=True):
                st.balloons()
                st.success("🎉 Analysis started! Check the tabs below for results.")
        
    # Main content area
    if uploaded_files:
        # Initialize session state for results
        if 'analysis_results' not in st.session_state:
            st.session_state.analysis_results = {}
        
        # Process files with enhanced progress feedback
        progress_container = st.container()
        with progress_container:
            progress_bar = st.progress(0)
            status_text = st.empty()
            
        with st.spinner("🔬 Processing uploaded files..."):
            status_text.text("📁 Reading and validating files...")
            progress_bar.progress(25)
            try:
                processed_files = parse_uploaded_files(uploaded_files)
                
                if processed_files:
                    progress_bar.progress(100)
                    status_text.text("✅ Processing complete!")
                    st.success(f"🎉 Successfully loaded {len(processed_files)} files with {sum(len(data['sequences']) for data in processed_files.values())} total sequences")
                    
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
    
    # Hero section
    st.markdown("""
    <div style="text-align: center; padding: 2rem; background: linear-gradient(135deg, #262730, #1e2130); border-radius: 1rem; margin: 2rem 0; border: 1px solid #00d4aa;">
        <h3 style="color: #00d4aa; margin-bottom: 1rem;">🧬 Professional Genome Analysis Tool</h3>
        <p style="color: #fafafa; font-size: 1.1rem; margin-bottom: 1.5rem;">Upload your DNA/RNA sequences and get comprehensive analysis with interactive visualizations</p>
        <p style="color: #00d4aa; font-weight: bold;">✨ Ready for deployment at arditmishra.com ✨</p>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        ### 🚀 **Advanced Features**
        ✅ **Multi-format Support**: FASTA and FASTQ files  
        ✅ **Real-time Analysis**: GC content, nucleotide composition, quality metrics  
        ✅ **Interactive Visualizations**: Dynamic charts with Plotly  
        ✅ **Batch Processing**: Handle up to 500 sequences efficiently  
        ✅ **Pattern Search**: Find specific DNA motifs and sequences  
        ✅ **Export Options**: CSV, JSON, and comprehensive text reports  
        ✅ **Dark Theme**: Easy on the eyes for long analysis sessions
        """)
    
    with col2:
        st.markdown("""
        ### 📋 **Quick Start Guide**
        
        **Step 1:** 📁 Upload Files  
        Use the sidebar to drag & drop your sequence files
        
        **Step 2:** ⚙️ Configure Analysis  
        Choose which analyses to run and set parameters
        
        **Step 3:** 🚀 Start Processing  
        Click "Start Analysis" and watch the magic happen
        
        **Step 4:** 📊 Explore Results  
        Browse through Overview, Analysis, and Visualizations
        
        **Step 5:** 📥 Download Reports  
        Export your results in multiple formats
        """)
    
    # Interactive demo section
    st.markdown('<h3 class="section-header">🧪 Ready to Analyze Your Data?</h3>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("""
        <div class="metric-card">
            <h4 style="color: #00d4aa; margin-bottom: 0.5rem;">📁 FASTA Files</h4>
            <p style="color: #fafafa; margin: 0;">DNA/RNA sequences with headers</p>
            <code style="color: #00d4aa;">>.fasta, .fa</code>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="metric-card">
            <h4 style="color: #00d4aa; margin-bottom: 0.5rem;">📊 FASTQ Files</h4>
            <p style="color: #fafafa; margin: 0;">Sequences with quality scores</p>
            <code style="color: #00d4aa;">.fastq, .fq</code>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown("""
        <div class="metric-card">
            <h4 style="color: #00d4aa; margin-bottom: 0.5rem;">🚀 Batch Mode</h4>
            <p style="color: #fafafa; margin: 0;">Process multiple files at once</p>
            <code style="color: #00d4aa;">Up to 50 files</code>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("""
    <div style="text-align: center; padding: 1.5rem; background: linear-gradient(135deg, #00d4aa20, #00b89420); border-radius: 1rem; margin: 2rem 0; border: 1px solid #00d4aa;">
        <p style="color: #fafafa; font-size: 1.1rem; margin: 0;">👈 <strong>Start by uploading your files in the sidebar</strong></p>
    </div>
    """, unsafe_allow_html=True)

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
    
    # Download options with enhanced UI
    st.markdown("### 📥 **Export Your Results**")
    st.markdown("Choose from multiple export formats to save your analysis results:")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        <div class="metric-card" style="text-align: center;">
            <h4 style="color: #00d4aa; margin-bottom: 0.5rem;">📊 CSV Format</h4>
            <p style="color: #fafafa; margin-bottom: 1rem;">Spreadsheet-ready data</p>
        </div>
        """, unsafe_allow_html=True)
        
        csv_data = generate_report(results, 'csv')
        st.download_button(
            label="📊 Download CSV Report",
            data=csv_data,
            file_name=f"genome_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv",
            use_container_width=True,
            type="secondary"
        )
    
    with col2:
        st.markdown("""
        <div class="metric-card" style="text-align: center;">
            <h4 style="color: #00d4aa; margin-bottom: 0.5rem;">📋 JSON Format</h4>
            <p style="color: #fafafa; margin-bottom: 1rem;">Structured data for APIs</p>
        </div>
        """, unsafe_allow_html=True)
        
        json_data = generate_report(results, 'json')
        st.download_button(
            label="📋 Download JSON Report",
            data=json_data,
            file_name=f"genome_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
            mime="application/json",
            use_container_width=True,
            type="secondary"
        )
    
    with col3:
        st.markdown("""
        <div class="metric-card" style="text-align: center;">
            <h4 style="color: #00d4aa; margin-bottom: 0.5rem;">📄 Full Report</h4>
            <p style="color: #fafafa; margin-bottom: 1rem;">Comprehensive text summary</p>
        </div>
        """, unsafe_allow_html=True)
        
        full_report = generate_report(results, 'full')
        st.download_button(
            label="📄 Download Full Report",
            data=full_report,
            file_name=f"genome_analysis_full_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
            mime="text/plain",
            use_container_width=True,
            type="secondary"
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
