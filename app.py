import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO, BytesIO
import zipfile
import json
from datetime import datetime
import hashlib
import time

from sequence_analyzer import SequenceAnalyzer
from visualizations import create_gc_content_plot, create_nucleotide_composition_plot, create_sequence_length_plot, create_quality_scores_plot, create_interactive_sequence_logo, create_comparative_analysis_plot
from utils import validate_file_format, parse_uploaded_files, generate_report, export_analysis_summary, create_batch_download_zip, validate_sequence_data

# Page configuration
st.set_page_config(
    page_title="Genome Sequencing Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Classic Scientific Dark Theme - Custom CSS
st.markdown("""
<style>
    /* Global Base */
    .stApp {
        background-color: #1A1D29;
        color: #E8EAED;
    }
    
    /* Main Header */
    .main-header {
        font-size: 2rem;
        font-weight: 600;
        color: #FFFFFF;
        margin: 0 0 0.5rem 0;
        padding: 0;
        letter-spacing: -0.5px;
    }
    
    .subtitle {
        font-size: 1rem;
        color: #9CA3AF;
        margin: 0 0 2rem 0;
        font-weight: 400;
    }
    
    /* Section Headers */
    .section-header {
        font-size: 1.25rem;
        font-weight: 600;
        color: #FFFFFF;
        margin: 2rem 0 1rem 0;
        padding-bottom: 0.5rem;
        border-bottom: 1px solid #373B4D;
    }
    
    /* Upload Zone - Prominent */
    .upload-zone {
        background: #25293C;
        border: 2px dashed #00C9A7;
        border-radius: 12px;
        padding: 3rem;
        text-align: center;
        margin: 2rem 0;
        transition: all 0.3s ease;
    }
    
    .upload-zone:hover {
        border-color: #00FFD1;
        background: #2A2E42;
    }
    
    .stFileUploader {
        background: transparent !important;
    }
    
    .stFileUploader > div {
        border: none !important;
        padding: 0 !important;
    }
    
    /* Cards */
    .metric-card {
        background: #25293C !important;
        padding: 1.5rem !important;
        border-radius: 8px !important;
        border: 1px solid #373B4D !important;
        margin-bottom: 1rem;
    }
    
    div[data-testid="stMetric"], div[data-testid="metric-container"] {
        background: #25293C !important;
        border-radius: 8px !important;
        border: 1px solid #373B4D !important;
        padding: 1rem !important;
    }
    
    /* Buttons */
    .stButton > button, button[kind="primary"] {
        background: #00C9A7 !important;
        color: #1A1D29 !important;
        border: none !important;
        border-radius: 6px !important;
        padding: 0.625rem 1.5rem !important;
        font-weight: 600 !important;
        transition: all 0.2s ease !important;
    }
    
    .stButton > button:hover, button[kind="primary"]:hover {
        background: #00FFD1 !important;
        transform: translateY(-1px) !important;
    }
    
    button[kind="secondary"] {
        background: transparent !important;
        color: #E8EAED !important;
        border: 1px solid #373B4D !important;
    }
    
    button[kind="secondary"]:hover {
        border-color: #00C9A7 !important;
        background: #25293C !important;
    }
    
    /* Sidebar */
    section[data-testid="stSidebar"] {
        background: #1F2230 !important;
        border-right: 1px solid #373B4D !important;
    }
    
    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        gap: 0;
        border-bottom: 1px solid #373B4D;
        background: transparent;
    }
    
    .stTabs [data-baseweb="tab"] {
        background: transparent;
        border: none;
        border-bottom: 2px solid transparent;
        padding: 0.75rem 1.5rem;
        color: #9CA3AF;
        font-weight: 500;
    }
    
    .stTabs [data-baseweb="tab"]:hover {
        color: #E8EAED;
        background: #25293C;
    }
    
    .stTabs [aria-selected="true"] {
        color: #00C9A7;
        border-bottom-color: #00C9A7;
        background: #25293C;
    }
    
    /* Progress Bar */
    .stProgress > div > div > div {
        background: #00C9A7 !important;
    }
    
    /* Expanders */
    .stExpander {
        background: #25293C;
        border: 1px solid #373B4D;
        border-radius: 8px;
        margin: 1rem 0;
    }
    
    /* Inputs */
    .stSelectbox > div > div, .stTextInput > div > div > input, .stNumberInput > div > div > input {
        background: #25293C !important;
        border: 1px solid #373B4D !important;
        border-radius: 6px !important;
        color: #E8EAED !important;
    }
    
    .stSelectbox > div > div:hover, .stTextInput > div > div > input:focus {
        border-color: #00C9A7 !important;
    }
    
    /* Checkboxes */
    .stCheckbox {
        color: #E8EAED !important;
    }
    
    /* Scrollbar */
    ::-webkit-scrollbar {
        width: 8px;
        height: 8px;
    }
    
    ::-webkit-scrollbar-track {
        background: #1A1D29;
    }
    
    ::-webkit-scrollbar-thumb {
        background: #373B4D;
        border-radius: 4px;
    }
    
    ::-webkit-scrollbar-thumb:hover {
        background: #4A4E63;
    }
    
    /* DNA Colors */
    .nucleotide-a { color: #FF6B9D; }
    .nucleotide-t { color: #FFC252; }
    .nucleotide-c { color: #00C9A7; }
    .nucleotide-g { color: #4D96FF; }
</style>
""", unsafe_allow_html=True)

def main():
    # Main header
    st.markdown('<h1 class="main-header">Genome Sequencing Analyzer</h1>', unsafe_allow_html=True)
    st.markdown('<p class="subtitle">Analyze DNA sequences for GC content, k-mers, and reading frames</p>', unsafe_allow_html=True)
    
    # File upload - always in main area
    st.markdown("---")
    
    # Upload section with prominent styling
    st.markdown("""
    <div style="background: #25293C; border: 2px dashed #00C9A7; border-radius: 12px; padding: 2.5rem; text-align: center; margin: 1.5rem 0;">
        <h3 style="color: #00C9A7; margin: 0 0 0.5rem 0; font-size: 1.5rem;">📁 Upload Your Sequence Files</h3>
        <p style="color: #9CA3AF; margin: 0 0 1.5rem 0; font-size: 1rem;">Drag and drop files here or click below to browse</p>
        <p style="color: #E8EAED; font-size: 0.9rem; margin: 0;">Supported: FASTA (.fasta, .fa) • FASTQ (.fastq, .fq)</p>
    </div>
    """, unsafe_allow_html=True)
    
    uploaded_files = st.file_uploader(
        "Choose files",
        type=['fasta', 'fa', 'fastq', 'fq'],
        accept_multiple_files=True,
        key="file_uploader",
        label_visibility="collapsed"
    )
    
    if not uploaded_files:
        # Information tabs
        st.markdown("---")
        
        tab1, tab2, tab3, tab4 = st.tabs(["📋 Overview", "🔬 Features & Capabilities", "📁 File Formats", "💡 Use Cases"])
        
        with tab1:
            st.markdown("### About This Tool")
            st.markdown("""
            The Genome Sequencing Analyzer is a comprehensive bioinformatics tool for analyzing DNA sequence data. 
            Upload your FASTA or FASTQ files to perform various analyses including composition, quality metrics, 
            and pattern recognition.
            """)
            
            st.markdown("### Quick Start")
            st.markdown("""
            1. **Upload** your sequence files using the upload area above
            2. **Select** the analyses you want to perform from the sidebar
            3. **View** interactive visualizations and detailed results
            4. **Export** your analysis reports for further use
            """)
            
            st.markdown("### Key Benefits")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("""
                - **Fast Processing** - Analyze sequences in seconds
                - **Interactive Visualizations** - Explore data with dynamic charts
                - **Batch Analysis** - Process multiple files simultaneously
                """)
            with col2:
                st.markdown("""
                - **Comprehensive Metrics** - GC content, k-mers, quality scores
                - **Export Ready** - Download results in various formats
                - **No Installation** - Web-based tool, works anywhere
                """)
        
        with tab2:
            st.markdown("### Core Analyses")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown("""
                <div class="metric-card">
                    <h4 style="color: #00C9A7; margin-bottom: 0.75rem;">🧬 GC Content Analysis</h4>
                    <p style="color: #9CA3AF; font-size: 0.9rem; line-height: 1.6;">
                    Calculate GC percentage for genome characterization. Essential for:
                    </p>
                    <ul style="color: #E8EAED; font-size: 0.85rem; margin-left: 1.2rem;">
                        <li>Species identification</li>
                        <li>Phylogenetic analysis</li>
                        <li>PCR primer design</li>
                        <li>Thermal stability prediction</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.markdown("""
                <div class="metric-card">
                    <h4 style="color: #00C9A7; margin-bottom: 0.75rem;">🔬 K-mer Frequency</h4>
                    <p style="color: #9CA3AF; font-size: 0.9rem; line-height: 1.6;">
                    Analyze k-mer patterns (3-8 nucleotides) for:
                    </p>
                    <ul style="color: #E8EAED; font-size: 0.85rem; margin-left: 1.2rem;">
                        <li>Sequence composition analysis</li>
                        <li>Motif discovery</li>
                        <li>Codon usage patterns</li>
                        <li>Repeat identification</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            with col3:
                st.markdown("""
                <div class="metric-card">
                    <h4 style="color: #00C9A7; margin-bottom: 0.75rem;">📖 Reading Frames</h4>
                    <p style="color: #9CA3AF; font-size: 0.9rem; line-height: 1.6;">
                    Six-frame translation analysis to identify:
                    </p>
                    <ul style="color: #E8EAED; font-size: 0.85rem; margin-left: 1.2rem;">
                        <li>Open reading frames (ORFs)</li>
                        <li>Potential coding regions</li>
                        <li>Start/stop codons</li>
                        <li>Protein translations</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            st.markdown("### Additional Analyses")
            st.markdown("""
            - **Nucleotide Composition** - A, T, C, G distribution with interactive charts
            - **Quality Scores** - Phred score analysis for FASTQ files  
            - **Pattern Search** - Find specific DNA motifs and sequences
            - **Batch Processing** - Analyze multiple files simultaneously
            - **Statistical Reports** - Comprehensive metrics and summaries
            """)
        
        with tab3:
            st.markdown("### Supported Formats")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("""
                #### FASTA Format (.fasta, .fa)
                
                **Description:**
                - Standard sequence format for DNA, RNA, or protein sequences
                - Header line starts with `>` followed by sequence identifier
                - Sequence data follows on subsequent lines
                
                **Example:**
                ```
                >sequence_1
                ATCGATCGATCG
                GCTAGCTAGCTA
                >sequence_2
                TTAACCGGTTAA
                ```
                
                **Use Cases:**
                - Reference genomes
                - Gene sequences
                - Protein sequences
                - Assembled contigs
                """)
            
            with col2:
                st.markdown("""
                #### FASTQ Format (.fastq, .fq)
                
                **Description:**
                - Sequencing data with quality scores
                - Four lines per sequence record
                - Includes Phred quality scores for each base
                
                **Example:**
                ```
                @SEQ_ID
                GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
                +
                !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
                ```
                
                **Use Cases:**
                - Illumina sequencing output
                - PacBio reads
                - Nanopore data
                - Quality control analysis
                """)
        
        with tab4:
            st.markdown("### Application Areas")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown("""
                #### 🔬 Research Applications
                
                **Genomics Research:**
                - Genome assembly validation
                - Gene prediction preprocessing
                - Comparative genomics studies
                - Metagenomics analysis
                
                **Molecular Biology:**
                - PCR primer design
                - Probe development
                - Restriction site analysis
                - Sequence optimization
                """)
            
            with col2:
                st.markdown("""
                #### 🏥 Clinical Applications
                
                **Diagnostics:**
                - Pathogen identification
                - Variant calling preparation
                - Quality control checks
                - Contamination detection
                
                **Precision Medicine:**
                - Biomarker discovery
                - Drug target analysis
                - Mutation screening
                - Resistance profiling
                """)
            
            with col3:
                st.markdown("""
                #### 📚 Educational Applications
                
                **Teaching:**
                - Bioinformatics courses
                - Sequence analysis labs
                - Demonstrating concepts
                - Hands-on training
                
                **Student Projects:**
                - Research projects
                - Data analysis assignments
                - Computational biology
                - Thesis work
                """)
        
        return
    
    st.success(f"✓ Loaded {len(uploaded_files)} file(s)")
    
    # Analysis options in sidebar
    with st.sidebar:
        st.markdown("### Analysis Options")
        
        analyze_gc_content = st.checkbox("GC Content", value=True)
        analyze_composition = st.checkbox("Nucleotide Composition", value=True)
        analyze_quality = st.checkbox("Quality Scores", value=True)
        
        st.markdown("---")
        st.markdown("### K-mer Analysis")
        analyze_kmers = st.checkbox("Enable K-mer Analysis", value=True)
        kmer_size = st.selectbox("K-mer size", [3, 4, 5, 6, 7, 8], index=0) if analyze_kmers else 3
        
        st.markdown("---")
        st.markdown("### Reading Frames")
        analyze_reading_frames = st.checkbox("Analyze Reading Frames", value=True)
        
        st.markdown("---")
        search_pattern = st.text_input("Search Pattern", placeholder="e.g., ATCG")
        
    # Main content area
    if uploaded_files:
        # Initialize session state for results and caching
        if 'analysis_results' not in st.session_state:
            st.session_state.analysis_results = {}
        if 'processed_files_cache' not in st.session_state:
            st.session_state.processed_files_cache = {}
        if 'analysis_params' not in st.session_state:
            st.session_state.analysis_params = {}
        
        # Process files with enhanced progress feedback
        progress_container = st.container()
        with progress_container:
            progress_bar = st.progress(0)
            status_text = st.empty()
            
        with st.spinner("Processing files..."):
            status_text.text("Reading files...")
            progress_bar.progress(25)
            try:
                processed_files = parse_uploaded_files(uploaded_files)
                
                if processed_files:
                    progress_bar.progress(100)
                    status_text.text("Complete")
                    st.success(f"Processed {len(processed_files)} files with {sum(len(data['sequences']) for data in processed_files.values())} sequences")
                    
                    
                    # Create tabs
                    tab1, tab2, tab3, tab4 = st.tabs([
                        "Overview", 
                        "Analysis", 
                        "Charts", 
                        "Export"
                    ])
                    
                    with tab1:
                        display_overview(processed_files)
                    
                    with tab2:
                        display_detailed_analysis(processed_files, analyze_gc_content, analyze_composition, analyze_quality, search_pattern, analyze_kmers, kmer_size, analyze_reading_frames)
                    
                    with tab3:
                        display_visualizations(processed_files, analyze_gc_content, analyze_composition, analyze_quality)
                    
                    with tab4:
                        display_download_section(processed_files)
                
                else:
                    st.error("No valid sequence files found. Please check your file formats.")
                    
            except Exception as e:
                st.error(f"Error processing files: {str(e)}")
                st.info("Please ensure your files are in valid FASTA or FASTQ format.")
    
def display_welcome_screen_legacy():
    """Display welcome screen with instructions"""
    st.markdown('<h1 class="main-header">Genome Sequencing Analyzer</h1>', unsafe_allow_html=True)
    
    # Simple hero section
    st.markdown("""
    <div style="text-align: center; padding: 2rem; background: #F8F9FA; border-radius: 8px; margin: 2rem 0; border: 1px solid #E0E0E0;">
        <p style="color: #2C3E50; font-size: 1.2rem; margin-bottom: 1rem;">Upload your DNA/RNA sequences for comprehensive genomic analysis</p>
        <p style="color: #7F8C8D; font-size: 1rem;">Analyze GC content, nucleotide composition, quality metrics, and more</p>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        ### Features
        - **Multi-format Support**: FASTA and FASTQ files  
        - **Real-time Analysis**: GC content, nucleotide composition, quality metrics  
        - **Interactive Charts**: Dynamic visualizations with Plotly  
        - **Batch Processing**: Handle multiple sequences efficiently  
        - **Pattern Search**: Find specific DNA motifs  
        - **Export Options**: CSV, JSON, and text reports
        """)
    
    with col2:
        st.markdown("""
        ### Quick Start
        
        **1. Upload Files**  
        Use the sidebar to upload your sequence files
        
        **2. Configure Analysis**  
        Choose which analyses to run
        
        **3. Start Processing**  
        Click "Start Analysis"
        
        **4. Explore Results**  
        Browse through the tabs below
        
        **5. Download Reports  
        Export your results in multiple formats
        """)
    
    # File format info
    st.markdown("### Supported Formats")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("""
        <div class="metric-card">
            <h4 style="color: #2C3E50; margin-bottom: 0.75rem;">📁 FASTA Files</h4>
            <p style="color: #7F8C8D; margin: 0 0 0.5rem 0;">DNA/RNA sequences with headers</p>
            <code style="color: #3498DB; background: #EBF5FB; padding: 0.25rem 0.5rem; border-radius: 4px;">.fasta, .fa</code>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="metric-card">
            <h4 style="color: #2C3E50; margin-bottom: 0.75rem;">📊 FASTQ Files</h4>
            <p style="color: #7F8C8D; margin: 0 0 0.5rem 0;">Sequences with quality scores</p>
            <code style="color: #3498DB; background: #EBF5FB; padding: 0.25rem 0.5rem; border-radius: 4px;">.fastq, .fq</code>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown("""
        <div class="metric-card">
            <h4 style="color: #2C3E50; margin-bottom: 0.75rem;">🚀 Batch Mode</h4>
            <p style="color: #7F8C8D; margin: 0 0 0.5rem 0;">Process multiple files at once</p>
            <code style="color: #3498DB; background: #EBF5FB; padding: 0.25rem 0.5rem; border-radius: 4px;">Up to 50 files</code>
        </div>
        """, unsafe_allow_html=True)
    
    st.info("👈 Start by uploading your files in the sidebar")

def display_overview(processed_files):
    """Display comprehensive overview with research insights"""
    st.markdown('<h2 class="section-header">Dataset Overview</h2>', unsafe_allow_html=True)
    
    if not processed_files:
        st.warning("No files to analyze")
        return
    
    # Calculate comprehensive statistics
    analyzer = SequenceAnalyzer()
    total_sequences = sum(len(data['sequences']) for data in processed_files.values())
    total_nucleotides = sum(
        sum(len(seq.seq) for seq in data['sequences']) 
        for data in processed_files.values()
    )
    file_types = [data['format'] for data in processed_files.values()]
    
    # Calculate GC content across all files
    all_gc_contents = []
    for data in processed_files.values():
        for seq in data['sequences']:
            all_gc_contents.append(analyzer.calculate_gc_content(str(seq.seq)))
    
    avg_gc_content = np.mean(all_gc_contents) if all_gc_contents else 0
    avg_length = total_nucleotides / total_sequences if total_sequences > 0 else 0
    
    # Enhanced metrics display
    st.markdown("### 📊 Dataset Statistics")
    col1, col2, col3, col4, col5 = st.columns(5)
    
    with col1:
        st.metric(
            "📁 Files Analyzed", 
            len(processed_files),
            help="Number of sequence files processed"
        )
    
    with col2:
        st.metric(
            "🧬 Total Sequences", 
            f"{total_sequences:,}",
            help="Total number of individual sequences across all files"
        )
    
    with col3:
        st.metric(
            "📏 Total Nucleotides", 
            f"{total_nucleotides:,} bp",
            help="Combined length of all sequences"
        )
    
    with col4:
        st.metric(
            "⚖️ Average GC Content", 
            f"{avg_gc_content:.1f}%",
            help="Overall GC content across the dataset"
        )
    
    with col5:
        st.metric(
            "📐 Average Length", 
            f"{avg_length:.0f} bp",
            help="Mean sequence length across all sequences"
        )
    
    # Research insights section
    st.markdown("### 🔬 Research Insights & Recommendations")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### 💡 Dataset Characteristics")
        
        # Genome type prediction
        if 30 <= avg_gc_content <= 45:
            st.info("🦠 **Mammalian-like GC content** - Suitable for vertebrate genome analysis")
        elif 45 <= avg_gc_content <= 70:
            st.info("🦠 **Bacterial-like GC content** - Good for prokaryotic studies")
        elif avg_gc_content > 70:
            st.info("🌱 **High-GC content** - May be plant DNA or actinobacteria")
        else:
            st.warning("🔬 **Unusual GC content** - Check for contamination or specialized organisms")
        
        # Dataset size assessment
        if total_nucleotides > 1000000:
            st.success("📊 **Large dataset** - Excellent for comprehensive analysis")
        elif total_nucleotides > 100000:
            st.info("📊 **Medium dataset** - Good for standard analysis")
        else:
            st.warning("📊 **Small dataset** - Limited analysis options")
    
    with col2:
        st.markdown("#### 🎯 Recommended Analysis")
        
        recommendations = []
        
        if total_sequences > 50:
            recommendations.append("🧬 Motif discovery and pattern analysis")
        if avg_length > 1000:
            recommendations.append("🔬 ORF prediction and gene analysis")
        if len(processed_files) > 1:
            recommendations.append("📊 Comparative genomics analysis")
        if total_nucleotides > 50000:
            recommendations.append("🌳 Phylogenetic and evolutionary analysis")
        
        if recommendations:
            for rec in recommendations:
                st.success(rec)
        else:
            st.info("💡 Upload more or longer sequences for advanced analysis")
    
    # Enhanced file details
    st.markdown("### 📄 Detailed File Analysis")
    
    for filename, data in processed_files.items():
        with st.expander(f"📁 {filename} - Analysis Summary", expanded=True):
            sequences = data['sequences']
            file_length = sum(len(seq.seq) for seq in sequences)
            file_avg_length = file_length / len(sequences) if sequences else 0
            
            # Calculate file-specific GC content
            file_gc_contents = [analyzer.calculate_gc_content(str(seq.seq)) for seq in sequences]
            file_avg_gc = np.mean(file_gc_contents) if file_gc_contents else 0
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Format", data['format'].upper())
                st.metric("Sequences", f"{len(sequences):,}")
            
            with col2:
                st.metric("Total Length", f"{file_length:,} bp")
                st.metric("Average Length", f"{file_avg_length:.1f} bp")
            
            with col3:
                st.metric("GC Content", f"{file_avg_gc:.1f}%")
                
                # Length category
                if file_avg_length > 5000:
                    length_cat = "🦣 Long sequences"
                elif file_avg_length > 1000:
                    length_cat = "🐘 Medium sequences"
                else:
                    length_cat = "🐭 Short sequences"
                st.write(length_cat)
            
            with col4:
                # File-specific recommendations
                file_recs = []
                if len(sequences) > 20:
                    file_recs.append("Statistical analysis")
                if file_avg_length > 2000:
                    file_recs.append("Structural analysis")
                if file_avg_gc > 60:
                    file_recs.append("High-GC organism study")
                
                if file_recs:
                    st.success(f"Best for: {', '.join(file_recs)}")
                else:
                    st.info("Basic analysis available")
            
            # Quick quality assessment
            quality_score = 100
            if len(sequences) < 10:
                quality_score -= 20
            if file_avg_length < 500:
                quality_score -= 15
            if file_avg_gc < 20 or file_avg_gc > 80:
                quality_score -= 10
            
            if quality_score >= 85:
                st.success(f"✅ Quality Score: {quality_score}/100")
            elif quality_score >= 65:
                st.warning(f"⚠️ Quality Score: {quality_score}/100")
            else:
                st.error(f"🚨 Quality Score: {quality_score}/100")
    
    # Data summary
    st.markdown("### 📈 Analysis Summary")
    st.success(f"✅ Dataset ready for comprehensive analysis with {len(processed_files)} files and {total_sequences:,} sequences")
    st.info("🔬 **Advanced features available**: Use the research analysis tabs for motif discovery, codon analysis, phylogenetic studies, and ORF prediction")

def display_detailed_analysis(processed_files, analyze_gc_content, analyze_composition, analyze_quality, search_pattern, analyze_kmers=True, kmer_size=3, analyze_reading_frames=True):
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
            
            # K-mer Analysis
            if analyze_kmers:
                st.markdown(f"**🧬 K-mer Analysis (k={kmer_size})**")
                try:
                    kmer_results = analyzer.analyze_kmers(sequences, k=kmer_size)
                    
                    # Display global k-mer statistics
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Total K-mers", f"{kmer_results['global_stats']['total_kmers']:,}")
                        st.metric("Unique K-mers", f"{kmer_results['global_stats']['unique_kmers']:,}")
                    
                    with col2:
                        # GC bias in k-mers
                        gc_bias = kmer_results['global_stats']['gc_bias']
                        st.metric("GC-rich K-mers", f"{gc_bias['gc_rich_percent']:.1f}%")
                        st.metric("AT-rich K-mers", f"{gc_bias['at_rich_percent']:.1f}%")
                    
                    with col3:
                        # Most common k-mers
                        most_common = kmer_results['global_stats']['most_common'][:5]
                        st.markdown("**Most Common:**")
                        for kmer, count in most_common:
                            st.write(f"{kmer}: {count:,}")
                    
                    # Show k-mer diversity
                    if kmer_results['per_sequence']:
                        avg_diversity = np.mean([seq['diversity'] for seq in kmer_results['per_sequence']])
                        st.info(f"📊 Average Shannon Diversity: {avg_diversity:.2f}")
                        
                        if avg_diversity > 3.0:
                            st.success("🌟 High k-mer diversity - good sequence complexity")
                        elif avg_diversity > 2.0:
                            st.info("📊 Moderate k-mer diversity")
                        else:
                            st.warning("⚠️ Low k-mer diversity - check for repeats or bias")
                    
                except Exception as e:
                    st.error(f"K-mer analysis failed: {str(e)}")
            
            # Reading Frame Analysis
            if analyze_reading_frames:
                st.markdown("**📖 Reading Frame Analysis**")
                try:
                    frame_results = analyzer.analyze_reading_frames(sequences)
                    
                    # Display reading frame statistics
                    for frame_data in frame_results['frame_stats']:
                        if frame_data['sequence_id'] == sequences[0].id:  # Show first sequence as example
                            st.markdown(f"**Analysis for: {frame_data['sequence_id']}**")
                            
                            # Create frame comparison table
                            frame_table = []
                            for frame_id, frame_info in frame_data['frames'].items():
                                frame_table.append({
                                    'Frame': frame_id,
                                    'Total Codons': frame_info['total_codons'],
                                    'Start Codons (ATG)': frame_info['start_codons'],
                                    'Stop Codons': frame_info['stop_codons'],
                                    'Longest ORF (bp)': frame_info['longest_orf'],
                                    'Stop Frequency': f"{frame_info['stop_frequency']:.3f}"
                                })
                            
                            frame_df = pd.DataFrame(frame_table)
                            st.dataframe(frame_df, use_container_width=True)
                            
                            # Find most promising reading frame
                            best_frame = max(frame_data['frames'].items(), 
                                           key=lambda x: x[1]['longest_orf'])
                            st.success(f"🎯 **Best reading frame**: {best_frame[0]} with longest ORF of {best_frame[1]['longest_orf']} bp")
                            
                            # Reading frame insights
                            total_frames = len(frame_data['frames'])
                            frames_with_orfs = sum(1 for f in frame_data['frames'].values() if f['longest_orf'] > 100)
                            
                            if frames_with_orfs > 3:
                                st.info("🧬 Multiple frames have long ORFs - potentially coding sequence")
                            elif frames_with_orfs > 0:
                                st.info("📊 Some reading frames show coding potential")
                            else:
                                st.warning("⚠️ No significant ORFs found - may be non-coding sequence")
                            
                            break  # Only show first sequence to avoid clutter
                    
                    # Summary for all sequences
                    if len(frame_results['frame_stats']) > 1:
                        all_longest_orfs = []
                        for frame_data in frame_results['frame_stats']:
                            max_orf = max(f['longest_orf'] for f in frame_data['frames'].values())
                            all_longest_orfs.append(max_orf)
                        
                        avg_longest_orf = np.mean(all_longest_orfs)
                        st.info(f"📈 Average longest ORF across all sequences: {avg_longest_orf:.1f} bp")
                    
                except Exception as e:
                    st.error(f"Reading frame analysis failed: {str(e)}")
            
            # Pattern search
            if search_pattern:
                st.markdown(f"**🔍 Pattern Search: '{search_pattern}'**")
                matches = analyzer.search_pattern(sequences, search_pattern)
                if matches:
                    match_df = pd.DataFrame(matches)
                    st.dataframe(match_df, hide_index=True)
                    
                    # Pattern analysis insights
                    st.info(f"🎯 Found {len(matches)} matches across {len(set(m['sequence_id'] for m in matches))} sequences")
                else:
                    st.info(f"Pattern '{search_pattern}' not found in any sequences")
                    st.markdown("**💡 Try these common patterns:**")
                    st.write("• ATG (start codon) • TAA, TAG, TGA (stop codons) • TATA (promoter) • CG (methylation sites)")

def display_visualizations(processed_files, analyze_gc_content, analyze_composition, analyze_quality):
    """Display interactive visualizations"""
    st.markdown('<h2 class="section-header">📈 Interactive Visualizations</h2>', unsafe_allow_html=True)
    
    # Add helpful guidance
    st.info("💡 **Tip**: All charts are interactive! Hover, zoom, and click to explore your data in detail.")
    
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
    
    # Add new advanced visualizations
    if len(processed_files) > 1:
        st.markdown("### 📊 File Comparison Dashboard")
        st.info("📋 Compare statistics across all uploaded files")
        comp_fig = create_comparative_analysis_plot(processed_files)
        st.plotly_chart(comp_fig, use_container_width=True)
    
    # Sequence logo for pattern analysis
    st.markdown("### 🧬 Sequence Logo Analysis")
    
    # File selector for sequence logo
    if len(processed_files) > 1:
        selected_file = st.selectbox(
            "Select file for sequence logo analysis",
            list(processed_files.keys()),
            help="Choose which file to analyze for nucleotide patterns"
        )
        sequences_for_logo = processed_files[selected_file]['sequences']
    else:
        selected_file = list(processed_files.keys())[0]
        sequences_for_logo = processed_files[selected_file]['sequences']
    
    # Sequence logo options
    col1, col2 = st.columns(2)
    with col1:
        max_sequences = st.slider(
            "Number of sequences to analyze", 
            1, 
            min(100, len(sequences_for_logo)), 
            min(50, len(sequences_for_logo)),
            help="More sequences = more accurate patterns, but slower processing"
        )
    with col2:
        if st.button("🔍 Generate Sequence Logo", type="secondary"):
            with st.spinner("Analyzing sequence patterns..."):
                logo_fig = create_interactive_sequence_logo(
                    sequences_for_logo[:max_sequences],
                    f"Sequence Logo - {selected_file}"
                )
                st.plotly_chart(logo_fig, use_container_width=True)
                st.success(f"✅ Analyzed {max_sequences} sequences from {selected_file}")

def display_download_section(processed_files):
    """Display enhanced download options"""
    st.markdown('<h2 class="section-header">📥 Download & Export Results</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    📊 **Export your analysis in multiple formats for different use cases:**
    - **CSV**: For spreadsheet analysis and data processing
    - **JSON**: For programmatic access and API integration  
    - **Text Report**: For documentation and sharing
    - **Complete Package**: All formats in one convenient ZIP file
    """)
    
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
    
    # Complete package download
    st.markdown("### 📦 **Complete Analysis Package**")
    st.markdown("🔥 **Most Popular!** Get everything in one ZIP file")
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("📦 Download Complete Package", type="primary", use_container_width=True):
            with st.spinner("Preparing your complete analysis package..."):
                zip_data = create_batch_download_zip(results)
                st.download_button(
                    label="📦 Download ZIP Package",
                    data=zip_data,
                    file_name=f"genome_analysis_complete_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
                    mime="application/zip",
                    use_container_width=True
                )
                st.success("✅ Package ready! Contains all reports, data, and documentation.")
    
    with col2:
        if st.button("📧 Generate Shareable Summary", use_container_width=True):
            summary = export_analysis_summary(results)
            st.text_area(
                "Copy and share this summary:",
                summary,
                height=200,
                help="Copy this text to share your results via email, chat, or documents"
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
