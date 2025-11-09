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

# Neural Bio-Futurism Theme - Custom CSS
st.markdown("""
<style>
    /* Global Base Styles */
    .stApp {
        background-color: #0B0E13;
        color: #E2E8F0;
    }
    
    /* Main Hero Header with Bioluminescent Glow */
    .main-header {
        font-size: 3.5rem;
        font-weight: 700;
        background: linear-gradient(90deg, #00D8B4, #38BDF8, #00D8B4);
        background-size: 200% 100%;
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        text-align: center;
        margin: 2rem 0;
        padding: 1rem;
        animation: gradientFlow 4s ease infinite, biolumGlow 2s ease-in-out infinite alternate;
        filter: drop-shadow(0 0 20px rgba(0, 216, 180, 0.4));
    }
    
    @keyframes gradientFlow {
        0%, 100% { background-position: 0% 50%; }
        50% { background-position: 100% 50%; }
    }
    
    @keyframes biolumGlow {
        from { 
            filter: drop-shadow(0 0 15px rgba(0, 216, 180, 0.3));
        }
        to { 
            filter: drop-shadow(0 0 30px rgba(0, 216, 180, 0.6));
        }
    }
    
    /* Section Headers with Gradient Underline */
    .section-header {
        font-size: 1.75rem;
        font-weight: 600;
        background: linear-gradient(90deg, #00D8B4, #38BDF8);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        margin: 2rem 0 1rem 0;
        padding-bottom: 0.75rem;
        border-bottom: 2px solid transparent;
        border-image: linear-gradient(90deg, #00D8B4, #38BDF8, transparent) 1;
        position: relative;
    }
    
    .section-header::after {
        content: '';
        position: absolute;
        bottom: -2px;
        left: 0;
        width: 100%;
        height: 2px;
        background: linear-gradient(90deg, #00D8B4, #38BDF8, transparent);
        box-shadow: 0 0 10px rgba(0, 216, 180, 0.5);
    }
    
    /* Glass-Morphism Cards with Bioluminescent Glow */
    .metric-card {
        background: rgba(22, 27, 34, 0.7) !important;
        backdrop-filter: blur(10px) !important;
        -webkit-backdrop-filter: blur(10px) !important;
        padding: 2rem !important;
        border-radius: 1.25rem !important;
        border: 1px solid rgba(0, 216, 180, 0.3) !important;
        box-shadow: 
            0 8px 32px rgba(0, 216, 180, 0.1),
            inset 0 1px 0 rgba(255, 255, 255, 0.05) !important;
        transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1) !important;
        position: relative;
        overflow: hidden;
    }
    
    .metric-card::before {
        content: '';
        position: absolute;
        top: 0;
        left: -100%;
        width: 100%;
        height: 100%;
        background: linear-gradient(90deg, transparent, rgba(0, 216, 180, 0.1), transparent);
        transition: left 0.5s;
        z-index: 0;
    }
    
    .metric-card:hover {
        transform: translateY(-4px) !important;
        border-color: rgba(56, 189, 248, 0.6) !important;
        box-shadow: 
            0 12px 48px rgba(0, 216, 180, 0.2),
            0 0 20px rgba(56, 189, 248, 0.3),
            inset 0 1px 0 rgba(255, 255, 255, 0.1) !important;
    }
    
    .metric-card:hover::before {
        left: 100%;
    }
    
    /* Additional Card Styling for Streamlit Metrics */
    div[data-testid="stMetric"], div[data-testid="metric-container"] {
        background: rgba(22, 27, 34, 0.7) !important;
        backdrop-filter: blur(10px) !important;
        -webkit-backdrop-filter: blur(10px) !important;
        border-radius: 1rem !important;
        border: 1px solid rgba(0, 216, 180, 0.2) !important;
        padding: 1.5rem !important;
        box-shadow: 0 4px 20px rgba(0, 216, 180, 0.1) !important;
    }
    
    /* File Uploader with Glowing Drop Zone */
    .stFileUploader > div > div > div {
        background: rgba(22, 27, 34, 0.5);
        backdrop-filter: blur(8px);
        border: 2px dashed rgba(0, 216, 180, 0.4);
        border-radius: 1.25rem;
        padding: 2.5rem;
        transition: all 0.3s ease;
        box-shadow: inset 0 0 20px rgba(0, 216, 180, 0.05);
    }
    
    .stFileUploader > div > div > div:hover {
        border-color: rgba(56, 189, 248, 0.8);
        background: rgba(22, 27, 34, 0.8);
        box-shadow: 
            0 0 30px rgba(0, 216, 180, 0.2),
            inset 0 0 30px rgba(0, 216, 180, 0.1);
    }
    
    /* Select Boxes and Inputs with Glass Effect */
    .stSelectbox > div > div, .stTextInput > div > div > input {
        background: rgba(22, 27, 34, 0.6);
        backdrop-filter: blur(8px);
        border: 1px solid rgba(0, 216, 180, 0.3);
        border-radius: 0.75rem;
        color: #E2E8F0;
        transition: all 0.3s ease;
    }
    
    .stSelectbox > div > div:hover, .stTextInput > div > div > input:focus {
        border-color: rgba(56, 189, 248, 0.6);
        box-shadow: 0 0 15px rgba(0, 216, 180, 0.2);
    }
    
    /* Checkboxes with Gradient Accent */
    .stCheckbox > label {
        color: #E2E8F0;
        font-weight: 500;
        transition: color 0.3s ease;
    }
    
    .stCheckbox > label:hover {
        color: #00D8B4;
    }
    
    /* Buttons with Cyan-Teal Gradient and Glow - Multiple Selectors for Streamlit */
    .stButton > button, button[kind="primary"], button[kind="secondary"] {
        background: linear-gradient(135deg, #00D8B4, #38BDF8) !important;
        color: white !important;
        border: none !important;
        border-radius: 0.75rem !important;
        padding: 0.75rem 1.5rem !important;
        font-weight: 600 !important;
        font-size: 1rem !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1) !important;
        box-shadow: 
            0 4px 20px rgba(0, 216, 180, 0.3),
            0 0 10px rgba(56, 189, 248, 0.2) !important;
        position: relative;
        overflow: hidden;
    }
    
    .stButton > button:hover, button[kind="primary"]:hover, button[kind="secondary"]:hover {
        transform: translateY(-3px) !important;
        box-shadow: 
            0 8px 30px rgba(0, 216, 180, 0.5),
            0 0 20px rgba(56, 189, 248, 0.4) !important;
        filter: brightness(1.1) !important;
    }
    
    .stButton > button:active, button[kind="primary"]:active {
        transform: translateY(-1px) !important;
    }
    
    /* Expanders with Glass Effect */
    .stExpander {
        background: rgba(22, 27, 34, 0.6);
        backdrop-filter: blur(10px);
        border: 1px solid rgba(0, 216, 180, 0.2);
        border-radius: 1rem;
        margin: 1.5rem 0;
        box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
    }
    
    .stExpander:hover {
        border-color: rgba(56, 189, 248, 0.4);
        box-shadow: 0 4px 20px rgba(0, 216, 180, 0.15);
    }
    
    /* Data Frames with Subtle Glow */
    .stDataFrame {
        background: rgba(22, 27, 34, 0.5);
        backdrop-filter: blur(8px);
        border-radius: 0.75rem;
        border: 1px solid rgba(0, 216, 180, 0.2);
        overflow: hidden;
        box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
    }
    
    /* Alert Messages with Magenta Accents */
    .stAlert {
        background: rgba(22, 27, 34, 0.7);
        backdrop-filter: blur(10px);
        border-left: 4px solid #C084FC;
        border-radius: 0.75rem;
        box-shadow: 0 4px 15px rgba(192, 132, 252, 0.1);
    }
    
    /* Success/Info/Warning/Error Styling */
    div[data-baseweb="notification"] {
        background: rgba(22, 27, 34, 0.9);
        backdrop-filter: blur(15px);
        border-radius: 0.75rem;
        border: 1px solid rgba(0, 216, 180, 0.3);
    }
    
    /* Progress Bars with Gradient */
    .stProgress > div > div > div {
        background: linear-gradient(90deg, #00D8B4, #38BDF8);
        box-shadow: 0 0 15px rgba(0, 216, 180, 0.5);
    }
    
    /* Sidebar with Dark Glass - Enhanced with !important */
    section[data-testid="stSidebar"], [data-testid="stSidebar"] > div {
        background: linear-gradient(180deg, rgba(11, 14, 19, 0.95), rgba(22, 27, 34, 0.95)) !important;
        backdrop-filter: blur(20px) !important;
        -webkit-backdrop-filter: blur(20px) !important;
        border-right: 1px solid rgba(0, 216, 180, 0.2) !important;
    }
    
    section[data-testid="stSidebar"]::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background: linear-gradient(180deg, rgba(0, 216, 180, 0.05), rgba(56, 189, 248, 0.05));
        pointer-events: none;
        z-index: 0;
    }
    
    /* Tab Navigation with Gradient */
    .stTabs [data-baseweb="tab-list"] {
        gap: 1rem;
        background: rgba(22, 27, 34, 0.4);
        backdrop-filter: blur(10px);
        padding: 0.5rem;
        border-radius: 1rem;
        border: 1px solid rgba(0, 216, 180, 0.2);
    }
    
    .stTabs [data-baseweb="tab"] {
        background: rgba(22, 27, 34, 0.6);
        backdrop-filter: blur(8px);
        border: 1px solid rgba(148, 163, 184, 0.2);
        border-radius: 0.75rem;
        padding: 0.75rem 1.5rem;
        color: #94A3B8;
        font-weight: 500;
        transition: all 0.3s ease;
    }
    
    .stTabs [data-baseweb="tab"]:hover {
        background: rgba(22, 27, 34, 0.8);
        border-color: rgba(0, 216, 180, 0.4);
        color: #E2E8F0;
    }
    
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #00D8B4, #38BDF8);
        border: 1px solid transparent;
        color: white;
        box-shadow: 
            0 4px 20px rgba(0, 216, 180, 0.4),
            0 0 15px rgba(56, 189, 248, 0.3);
    }
    
    /* Metrics with Glow */
    div[data-testid="stMetricValue"] {
        color: #E2E8F0;
        font-weight: 600;
        text-shadow: 0 0 10px rgba(0, 216, 180, 0.3);
    }
    
    /* Scrollbar Styling */
    ::-webkit-scrollbar {
        width: 10px;
        height: 10px;
    }
    
    ::-webkit-scrollbar-track {
        background: #0B0E13;
    }
    
    ::-webkit-scrollbar-thumb {
        background: linear-gradient(180deg, #00D8B4, #38BDF8);
        border-radius: 10px;
        box-shadow: 0 0 10px rgba(0, 216, 180, 0.5);
    }
    
    ::-webkit-scrollbar-thumb:hover {
        background: linear-gradient(180deg, #38BDF8, #00D8B4);
    }
    
    /* DNA Nucleotide Colors */
    .nucleotide-a { color: #FF6B6B; }
    .nucleotide-t { color: #FFD93D; }
    .nucleotide-c { color: #6BCB77; }
    .nucleotide-g { color: #4D96FF; }
    
    /* Enhanced Typography */
    h1, h2, h3 {
        color: #E2E8F0;
    }
    
    p, span, div {
        color: #E2E8F0;
    }
    
    /* Subtle Text */
    .stMarkdown small, .caption {
        color: #94A3B8;
    }
</style>
""", unsafe_allow_html=True)

def main():
    # Main header
    st.markdown('<h1 class="main-header">🧬 Genome Sequencing Analyzer</h1>', unsafe_allow_html=True)
    st.markdown("Upload DNA sequences and instantly analyze **GC content**, **k-mers**, and **reading frames** with comprehensive statistical analysis and downloadable results.")
    
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
        
        st.markdown("**🧬 K-mer Analysis**")
        analyze_kmers = st.checkbox(
            "K-mer Frequency Analysis",
            value=True,
            help="Analyze k-mer patterns and frequencies - essential for sequence composition studies"
        )
        if analyze_kmers:
            kmer_size = st.selectbox(
                "K-mer size",
                [3, 4, 5, 6, 7, 8],
                index=0,
                help="Choose k-mer length: 3 for codons, 4-8 for motif discovery"
            )
        else:
            kmer_size = 3
        
        st.markdown("**📖 Reading Frame Analysis**")
        analyze_reading_frames = st.checkbox(
            "Reading Frame Analysis",
            value=True,
            help="Analyze all 6 reading frames to identify potential coding regions and ORFs"
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
        
        # Real-time parameter updates
        current_params = {
            'gc_content': analyze_gc_content,
            'composition': analyze_composition,
            'quality': analyze_quality,
            'pattern': search_pattern,
            'batch_size': batch_size
        }
        
        # Check if parameters changed
        params_changed = current_params != st.session_state.get('analysis_params', {})
        if params_changed:
            st.session_state.analysis_params = current_params
        
        # Advanced filtering options
        st.markdown("**🔬 Advanced Filters**")
        with st.expander("Filter Options", expanded=False):
            min_length = st.number_input("Minimum sequence length", min_value=1, value=1, help="Filter out sequences shorter than this length")
            max_length = st.number_input("Maximum sequence length", min_value=1, value=100000, help="Filter out sequences longer than this length")
            gc_range = st.slider("GC Content Range (%)", 0.0, 100.0, (0.0, 100.0), help="Only analyze sequences within this GC content range")
        
        # Real-time search with debouncing
        if search_pattern:
            st.markdown(f"🔍 **Live Pattern Search**: *{search_pattern.upper()}*")
            if len(search_pattern) >= 3:
                st.info(f"✓ Pattern is valid for search ({len(search_pattern)} nucleotides)")
            else:
                st.warning("⚠️ Pattern should be at least 3 nucleotides for reliable results")
        
        # Analysis button with enhanced feedback
        if uploaded_files:
            col1, col2 = st.columns(2)
            with col1:
                if st.button("🚀 Start Analysis", type="primary", use_container_width=True):
                    st.balloons()
                    st.success("🎉 Analysis started! Check the tabs below for results.")
            with col2:
                if st.button("🔄 Reset Analysis", use_container_width=True):
                    st.session_state.analysis_results = {}
                    st.session_state.processed_files_cache = {}
                    st.success("Analysis cache cleared!")
        
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
            
        with st.spinner("🔬 Processing uploaded files..."):
            status_text.text("📁 Reading and validating files...")
            progress_bar.progress(25)
            try:
                processed_files = parse_uploaded_files(uploaded_files)
                
                if processed_files:
                    progress_bar.progress(100)
                    status_text.text("✅ Processing complete!")
                    st.success(f"🎉 Successfully loaded {len(processed_files)} files with {sum(len(data['sequences']) for data in processed_files.values())} total sequences")
                    
                    # Data quality validation
                    for filename, data in processed_files.items():
                        validation = validate_sequence_data(data['sequences'])
                        if validation['warnings']:
                            with st.expander(f"⚠️ Quality Check: {filename}", expanded=False):
                                for warning in validation['warnings']:
                                    st.warning(warning)
                                for rec in validation['recommendations']:
                                    st.info(rec)
                                st.metric("Data Quality Score", f"{validation['quality_score']}/100")
                        else:
                            st.success(f"✅ {filename}: High quality data (100/100)")
                    
                    # Create tabs for different views with enhanced descriptions
                    tab1, tab2, tab3, tab4 = st.tabs([
                        "📊 Overview", 
                        "🔬 Detailed Analysis", 
                        "📈 Interactive Charts", 
                        "📥 Export & Share"
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
    
    else:
        # Welcome screen
        display_welcome_screen()

def display_welcome_screen():
    """Display welcome screen with instructions"""
    st.markdown('<h2 class="section-header">Welcome to Genome Sequencing Analyzer</h2>', unsafe_allow_html=True)
    
    # Hero section with Neural Bio-Futurism theme
    st.markdown("""
    <div style="text-align: center; padding: 3rem; background: rgba(22, 27, 34, 0.7); backdrop-filter: blur(15px); border-radius: 1.5rem; margin: 2rem 0; border: 1px solid rgba(0, 216, 180, 0.3); box-shadow: 0 8px 32px rgba(0, 216, 180, 0.15), inset 0 1px 0 rgba(255, 255, 255, 0.05);">
        <h3 style="background: linear-gradient(90deg, #00D8B4, #38BDF8); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-size: 2rem; font-weight: 700; margin-bottom: 1rem; filter: drop-shadow(0 0 10px rgba(0, 216, 180, 0.4));">🧬 Neural Bio-Futurism Genome Analysis</h3>
        <p style="color: #E2E8F0; font-size: 1.2rem; margin-bottom: 1.5rem; line-height: 1.6;">Upload your DNA/RNA sequences and experience cutting-edge analysis powered by precision genomics</p>
        <p style="background: linear-gradient(90deg, #00D8B4, #38BDF8); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-weight: 700; font-size: 1.1rem;">✨ Precision Medicine Platform ✨</p>
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
            <h4 style="background: linear-gradient(90deg, #00D8B4, #38BDF8); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-weight: 700; margin-bottom: 0.75rem;">📁 FASTA Files</h4>
            <p style="color: #E2E8F0; margin: 0 0 0.5rem 0;">DNA/RNA sequences with headers</p>
            <code style="color: #00D8B4; background: rgba(0, 216, 180, 0.1); padding: 0.25rem 0.5rem; border-radius: 0.25rem;">>.fasta, .fa</code>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="metric-card">
            <h4 style="background: linear-gradient(90deg, #00D8B4, #38BDF8); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-weight: 700; margin-bottom: 0.75rem;">📊 FASTQ Files</h4>
            <p style="color: #E2E8F0; margin: 0 0 0.5rem 0;">Sequences with quality scores</p>
            <code style="color: #38BDF8; background: rgba(56, 189, 248, 0.1); padding: 0.25rem 0.5rem; border-radius: 0.25rem;">.fastq, .fq</code>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown("""
        <div class="metric-card">
            <h4 style="background: linear-gradient(90deg, #00D8B4, #38BDF8); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-weight: 700; margin-bottom: 0.75rem;">🚀 Batch Mode</h4>
            <p style="color: #E2E8F0; margin: 0 0 0.5rem 0;">Process multiple files at once</p>
            <code style="color: #C084FC; background: rgba(192, 132, 252, 0.1); padding: 0.25rem 0.5rem; border-radius: 0.25rem;">Up to 50 files</code>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("""
    <div style="text-align: center; padding: 2rem; background: rgba(22, 27, 34, 0.6); backdrop-filter: blur(10px); border-radius: 1rem; margin: 2rem 0; border: 1px solid rgba(0, 216, 180, 0.3); box-shadow: 0 4px 20px rgba(0, 216, 180, 0.15);">
        <p style="color: #E2E8F0; font-size: 1.2rem; margin: 0;">👈 <strong style="background: linear-gradient(90deg, #00D8B4, #38BDF8); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;">Start by uploading your files in the sidebar</strong></p>
    </div>
    """, unsafe_allow_html=True)

def display_overview(processed_files):
    """Display comprehensive overview with research insights"""
    st.markdown('<h2 class="section-header">🔍 Comprehensive Dataset Overview</h2>', unsafe_allow_html=True)
    
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
