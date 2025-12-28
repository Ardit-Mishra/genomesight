"""
Genome Sequencing Analyzer

A comprehensive bioinformatics web application for DNA sequence analysis.
Provides GC content calculation, nucleotide composition, k-mer analysis,
ORF detection, motif discovery, and interactive visualizations.

Author: Ardit Mishra
Version: 1.0.0
"""

import streamlit as st
from io import StringIO
from datetime import datetime

from app.core.sequence_analyzer import SequenceAnalyzer
from app.core.io_fasta import parse_sequences, detect_format
from app.core.validation import validate_sequences
from app.core.orf import find_orfs, get_orf_summary
from app.core.motifs import search_motif, get_enzyme_list
from app.core.plots import (
    create_gc_plot,
    create_gc_sliding_window,
    create_composition_plot,
    create_kmer_plot,
    create_quality_plot,
    create_length_distribution,
    create_orf_plot
)
from app.core.export import generate_json_report, generate_csv_report, create_download_zip
from app.ui.styles import apply_dark_theme, COLORS
from app.ui.components import (
    display_sequence_card,
    display_metrics_row,
    display_analysis_results,
    file_uploader_section,
    display_orf_table,
    display_motif_results,
    display_stats_badges
)


st.set_page_config(
    page_title="Genome Sequencing Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


apply_dark_theme(st)


SAMPLE_FASTA = """>Human_BRCA1_partial
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGT
CCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATG
CTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGC
CTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGAC
>E_coli_16S_partial
AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAA
CAGGTCTTCGGACGCTGACGAGTGGCGAACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACT
ACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCG
>Synthetic_high_GC
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
CGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
"""


def initialize_session():
    """Initialize session state variables."""
    if 'analyzer' not in st.session_state:
        st.session_state.analyzer = SequenceAnalyzer()
    if 'sequences' not in st.session_state:
        st.session_state.sequences = None
    if 'file_format' not in st.session_state:
        st.session_state.file_format = None
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = None


class SampleFile:
    """Wrapper to mimic uploaded file interface for sample data."""
    
    def __init__(self, content: str, filename: str):
        self.name = filename
        self._content = content
    
    def getvalue(self):
        return self._content.encode('utf-8')
    
    def read(self):
        return self._content


def render_header():
    """Render the application header."""
    st.markdown("""
    <div style="text-align: center; padding: 2rem 0 1rem 0;">
        <h1 class="main-header">🧬 Genome Sequencing Analyzer</h1>
        <p class="subtitle">Comprehensive DNA/RNA sequence analysis platform</p>
    </div>
    """, unsafe_allow_html=True)
    
    display_stats_badges()


def render_info_tabs():
    """Render the information tabs section."""
    overview_tab, features_tab, formats_tab, cases_tab = st.tabs([
        "Overview", "Features & Capabilities", "File Formats", "Use Cases"
    ])
    
    with overview_tab:
        st.markdown("""
        This application provides a comprehensive suite of tools for analyzing 
        DNA and RNA sequences. Upload your sequence files in FASTA or FASTQ format 
        to get started with analysis.
        
        **Quick Start:**
        1. Upload your sequence file using the upload area below
        2. Or click "Try Sample File" to explore with example data
        3. View basic statistics and interactive visualizations
        4. Use the sidebar to access advanced analysis tools
        """)
    
    with features_tab:
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            **Basic Analysis**
            - GC content calculation
            - Nucleotide composition
            - Sequence length statistics
            - Quality score analysis (FASTQ)
            
            **Advanced Analysis**
            - K-mer frequency analysis
            - Open Reading Frame detection
            - Motif and pattern search
            - Restriction site mapping
            """)
        
        with col2:
            st.markdown("""
            **Visualization**
            - Interactive GC content plots
            - Sliding window analysis
            - Composition bar charts
            - Quality score distributions
            
            **Export Options**
            - JSON structured data
            - CSV spreadsheet format
            - Text reports
            - ZIP archives
            """)
    
    with formats_tab:
        st.markdown("""
        **FASTA Format**
        - Extension: `.fasta`, `.fa`, `.fna`
        - Contains sequence ID and nucleotide data
        - Widely used for reference genomes
        
        **FASTQ Format**
        - Extension: `.fastq`, `.fq`
        - Contains sequence + quality scores
        - Standard for sequencing reads
        
        **Example FASTA:**
        ```
        >sequence_id description
        ATGCATGCATGCATGC
        ```
        """)
    
    with cases_tab:
        st.markdown("""
        **Research Applications**
        - Genome annotation projects
        - Comparative genomics studies
        - Primer design workflows
        
        **Clinical Applications**
        - Variant analysis preparation
        - Quality control of sequencing runs
        - Pathogen identification
        
        **Educational Use**
        - Bioinformatics training
        - Sequence analysis demonstrations
        - Understanding genome structure
        """)


def render_upload_section():
    """Render the file upload section."""
    file_uploader_section()
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        uploaded_file = st.file_uploader(
            "Choose a sequence file",
            type=['fasta', 'fa', 'fna', 'fastq', 'fq'],
            label_visibility="collapsed"
        )
    
    with col2:
        if st.button("Try Sample File", type="secondary", use_container_width=True):
            st.session_state.sample_loaded = True
            uploaded_file = SampleFile(SAMPLE_FASTA, "sample_sequences.fasta")
    
    return uploaded_file


def process_uploaded_file(uploaded_file):
    """Process and analyze uploaded sequence file."""
    try:
        if hasattr(uploaded_file, 'getvalue'):
            content = uploaded_file.getvalue()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
        else:
            content = uploaded_file.read()
        
        filename = getattr(uploaded_file, 'name', 'uploaded_file')
        
        sequences, file_format = parse_sequences(content, filename)
        
        if not sequences:
            st.error("No valid sequences found in the file")
            return None
        
        st.session_state.sequences = sequences
        st.session_state.file_format = file_format
        
        analyzer = st.session_state.analyzer
        stats = analyzer.calculate_sequence_statistics(sequences)
        
        st.session_state.analysis_results = {
            'filename': filename,
            'format': file_format,
            **stats
        }
        
        return st.session_state.analysis_results
        
    except Exception as e:
        st.error(f"Error processing file: {str(e)}")
        return None


def render_basic_stats(results: dict):
    """Render basic statistics section."""
    st.markdown("### Analysis Results")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Sequences", f"{results.get('sequence_count', 0):,}")
    with col2:
        st.metric("Total Length", f"{results.get('total_length', 0):,} bp")
    with col3:
        st.metric("Avg GC%", f"{results.get('average_gc_content', 0):.1f}%")
    with col4:
        st.metric("Avg Length", f"{results.get('average_length', 0):.0f} bp")
    
    st.markdown("---")
    
    col1, col2 = st.columns(2)
    
    with col1:
        composition = results.get('nucleotide_composition', {})
        if composition:
            fig = create_composition_plot(composition)
            st.plotly_chart(fig, use_container_width=True)
    
    with col2:
        lengths = [len(seq.seq) for seq in st.session_state.sequences]
        if lengths:
            fig = create_length_distribution(lengths)
            st.plotly_chart(fig, use_container_width=True)


def render_sidebar():
    """Render the sidebar with analysis options."""
    with st.sidebar:
        st.markdown("### Analysis Options")
        
        st.markdown("---")
        
        st.markdown("#### K-mer Analysis")
        kmer_size = st.slider("K-mer size", min_value=2, max_value=10, value=3)
        run_kmer = st.button("Run K-mer Analysis")
        
        st.markdown("---")
        
        st.markdown("#### ORF Detection")
        min_orf_length = st.number_input("Minimum ORF length (bp)", 
                                          min_value=30, max_value=1000, value=100)
        include_reverse = st.checkbox("Include reverse complement", value=True)
        run_orf = st.button("Find ORFs")
        
        st.markdown("---")
        
        st.markdown("#### Motif Search")
        motif_pattern = st.text_input("Pattern (IUPAC)", placeholder="ATGR")
        enzyme = st.selectbox("Or select enzyme", 
                              ["Custom"] + list(get_enzyme_list().keys()))
        run_motif = st.button("Search Motif")
        
        st.markdown("---")
        
        st.markdown("#### Export")
        export_format = st.selectbox("Format", ["CSV", "JSON", "ZIP"])
        run_export = st.button("Export Results")
        
        return {
            'kmer': {'run': run_kmer, 'size': kmer_size},
            'orf': {'run': run_orf, 'min_length': min_orf_length, 
                    'reverse': include_reverse},
            'motif': {'run': run_motif, 'pattern': motif_pattern, 
                      'enzyme': enzyme},
            'export': {'run': run_export, 'format': export_format}
        }


def run_kmer_analysis(size: int):
    """Run k-mer analysis and display results."""
    sequences = st.session_state.sequences
    if not sequences:
        return
    
    with st.spinner("Analyzing k-mers..."):
        analyzer = st.session_state.analyzer
        kmer_results = analyzer.analyze_kmers(sequences, k=size)
        
        st.markdown(f"### {size}-mer Analysis")
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Total K-mers", f"{kmer_results['total_kmers']:,}")
        with col2:
            st.metric("Unique K-mers", f"{kmer_results['unique_kmers']:,}")
        
        fig = create_kmer_plot(kmer_results['all_kmers'], size)
        st.plotly_chart(fig, use_container_width=True)


def run_orf_analysis(min_length: int, include_reverse: bool):
    """Run ORF detection and display results."""
    sequences = st.session_state.sequences
    if not sequences:
        return
    
    with st.spinner("Finding ORFs..."):
        all_orfs = []
        for seq in sequences[:10]:
            orfs = find_orfs(str(seq.seq), min_length=min_length, 
                            include_reverse=include_reverse)
            all_orfs.extend(orfs)
        
        st.markdown("### Open Reading Frames")
        
        if all_orfs:
            summary = get_orf_summary(all_orfs)
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("ORFs Found", summary['total'])
            with col2:
                st.metric("Avg Length", f"{summary['average_length']:.0f} bp")
            with col3:
                st.metric("Longest", f"{summary['max_length']} bp")
            
            fig = create_orf_plot(all_orfs)
            st.plotly_chart(fig, use_container_width=True)
            
            display_orf_table(all_orfs)
        else:
            st.info("No ORFs found with the specified parameters")


def run_motif_analysis(pattern: str, enzyme: str):
    """Run motif search and display results."""
    sequences = st.session_state.sequences
    if not sequences:
        return
    
    if enzyme != "Custom":
        enzymes = get_enzyme_list()
        pattern = enzymes.get(enzyme, pattern)
    
    if not pattern:
        st.warning("Please enter a pattern or select an enzyme")
        return
    
    with st.spinner("Searching for motif..."):
        matches = search_motif(sequences, pattern)
        
        st.markdown(f"### Motif Search: {pattern}")
        display_motif_results(matches)


def run_export(format_type: str):
    """Export results in specified format."""
    results = st.session_state.analysis_results
    if not results:
        st.warning("No results to export")
        return
    
    filename = results.get('filename', 'sequences')
    data = {filename: results}
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    if format_type == "CSV":
        content = generate_csv_report(data)
        st.download_button(
            "Download CSV",
            content,
            f"analysis_{timestamp}.csv",
            "text/csv"
        )
    elif format_type == "JSON":
        content = generate_json_report(data)
        st.download_button(
            "Download JSON",
            content,
            f"analysis_{timestamp}.json",
            "application/json"
        )
    else:
        content = create_download_zip(data)
        st.download_button(
            "Download ZIP",
            content,
            f"analysis_{timestamp}.zip",
            "application/zip"
        )


def render_footer():
    """Render the application footer."""
    st.markdown("---")
    st.markdown(f"""
    <div style="text-align: center; padding: 1rem 0; color: {COLORS['text_secondary']};">
        <p style="margin: 0;">
            Genome Sequencing Analyzer v1.0.0 | Built with Streamlit
        </p>
        <p style="margin: 0.25rem 0 0 0; font-size: 0.85rem;">
            Ardit Mishra
        </p>
    </div>
    """, unsafe_allow_html=True)


def main():
    """Main application entry point."""
    initialize_session()
    
    render_header()
    
    if st.session_state.sequences is None:
        render_info_tabs()
    
    uploaded_file = render_upload_section()
    
    if uploaded_file is not None:
        results = process_uploaded_file(uploaded_file)
        
        if results:
            render_basic_stats(results)
            
            sidebar_options = render_sidebar()
            
            if sidebar_options['kmer']['run']:
                run_kmer_analysis(sidebar_options['kmer']['size'])
            
            if sidebar_options['orf']['run']:
                run_orf_analysis(
                    sidebar_options['orf']['min_length'],
                    sidebar_options['orf']['reverse']
                )
            
            if sidebar_options['motif']['run']:
                run_motif_analysis(
                    sidebar_options['motif']['pattern'],
                    sidebar_options['motif']['enzyme']
                )
            
            if sidebar_options['export']['run']:
                run_export(sidebar_options['export']['format'])
    
    render_footer()


if __name__ == "__main__":
    main()
