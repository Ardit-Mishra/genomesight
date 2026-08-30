"""
GenomeSight

A bioinformatics web application for DNA and RNA sequence analysis.
Provides GC content calculation, nucleotide composition, k-mer analysis,
ORF detection, motif discovery, pairwise alignment, and interactive
visualizations.

Author: Ardit Mishra
Version: 1.0.0
"""

import streamlit as st

# Must be the first Streamlit command the app executes. Keeping it directly
# below the import — rather than after the `app.*` imports — means adding a
# module-level st.* call to any of those modules can't break startup with
# "set_page_config() can only be called once, and must be called as the first
# Streamlit command".
st.set_page_config(
    page_title="GenomeSight",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

from io import StringIO
from datetime import datetime

from app.core.sequence_analyzer import SequenceAnalyzer
from app.core.io_fasta import parse_sequences, detect_format, extract_genbank_features
from app.core.validation import validate_sequences
from app.core.orf import find_orfs, get_orf_summary
from app.core.motifs import search_motif, get_enzyme_list
from app.core.alignment import align_sequences, format_alignment_display
from app.core.codons import count_codons, codon_usage_for_orf, merge_codon_usage, get_codon_table_rows
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
from app.ui.icons import icon, section_header
from app.ui.components import (
    display_sequence_card,
    display_metrics_row,
    display_analysis_results,
    file_uploader_section,
    render_empty_state,
    sidebar_section_header,
    display_orf_table,
    display_motif_results,
    display_stats_badges
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
    if 'sample_loaded' not in st.session_state:
        st.session_state.sample_loaded = False


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
    # Left-aligned: the sequence and its numbers are the subject, and a centred
    # masthead pushes them below the fold for no gain.
    st.markdown(f"""
    <div style="padding: 0.5rem 0 1.25rem 0; border-bottom: 1px solid var(--border, #232C36); margin-bottom: 1.25rem;">
        <div style="display: flex; align-items: center; gap: 0.6rem;">
            {icon('helix', 26, COLORS['text_primary'])}
            <h1 class="main-header" style="margin: 0;">GenomeSight</h1>
        </div>
        <p class="subtitle">DNA and RNA sequence analysis &middot; deterministic, no model</p>
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
            st.markdown("**Basic Analysis**")
            for name, text in (
                ("composition", "GC content & nucleotide composition"),
                ("ruler", "Sequence length statistics"),
            ):
                st.markdown(f'<div style="margin:.35rem 0 .35rem 0">{icon(name, 15, COLORS["text_secondary"])} &nbsp;{text}</div>', unsafe_allow_html=True)
            st.markdown('<div style="margin:.35rem 0 .35rem 0">Quality score analysis (FASTQ)</div>', unsafe_allow_html=True)

            st.markdown("**Advanced Analysis**")
            for name, text in (
                ("kmer", "K-mer frequency analysis"),
                ("readingFrame", "Open reading frame (ORF) detection"),
                ("readingFrame", "Codon usage tables — per-ORF and whole-sequence, with RSCU"),
                ("motif", "Motif & restriction-site search (IUPAC patterns)"),
                ("alignment", "Pairwise alignment — global (Needleman-Wunsch) or local (Smith-Waterman)"),
            ):
                st.markdown(f'<div style="margin:.35rem 0 .35rem 0">{icon(name, 15, COLORS["text_secondary"])} &nbsp;{text}</div>', unsafe_allow_html=True)

        with col2:
            st.markdown("""
            **Visualization**
            - GC content distribution across sequences
            - Nucleotide composition bar chart
            - Sequence length distribution
            - Quality score distribution (FASTQ only)

            **Export Options**
            - JSON structured data
            - CSV spreadsheet format
            - Text reports
            - ZIP archives
            """)
    
    with formats_tab:
        st.markdown(
            f'<div style="display:flex;align-items:center;gap:.5rem;margin-bottom:.4rem">'
            f'{icon("sequenceFile", 16, COLORS["text_primary"])}<strong>FASTA Format</strong></div>',
            unsafe_allow_html=True,
        )
        st.markdown("""
        - Extension: `.fasta`, `.fa`, `.fna`
        - Contains sequence ID and nucleotide data
        - Widely used for reference genomes
        """)
        st.markdown(
            f'<div style="display:flex;align-items:center;gap:.5rem;margin-bottom:.4rem">'
            f'{icon("sequenceFile", 16, COLORS["text_primary"])}<strong>FASTQ Format</strong></div>',
            unsafe_allow_html=True,
        )
        st.markdown("""
        - Extension: `.fastq`, `.fq`
        - Contains sequence + quality scores
        - Standard for sequencing reads

        **Example FASTA:**
        ```
        >sequence_id description
        ATGCATGCATGCATGC
        ```
        """)
        st.markdown(
            f'<div style="display:flex;align-items:center;gap:.5rem;margin-bottom:.4rem">'
            f'{icon("sequenceFile", 16, COLORS["text_primary"])}<strong>GenBank Format</strong></div>',
            unsafe_allow_html=True,
        )
        st.markdown("""
        - Extension: `.gb`, `.gbk`, `.genbank`
        - Sequence plus structured feature annotations (source, CDS, gene,
          etc.) — parsed via Biopython, features preserved
        - A file can hold multiple `LOCUS` records back to back
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
    if st.session_state.sequences is None:
        render_empty_state()

    file_uploader_section()

    col1, col2 = st.columns([3, 1])
    
    with col1:
        uploaded_file = st.file_uploader(
            "Choose a sequence file",
            type=['fasta', 'fa', 'fna', 'fastq', 'fq', 'gb', 'gbk', 'genbank'],
            label_visibility="collapsed"
        )
    
    with col2:
        if st.button("Try Sample File", type="secondary", use_container_width=True):
            st.session_state.sample_loaded = True

    # A real upload is retained by the file_uploader widget across reruns, but a
    # button is only True on the run that handled the click. Rebuilding the
    # sample from session state keeps it loaded — otherwise every sidebar action
    # after "Try Sample File" (k-mer, ORF, motif, alignment, export) reran with
    # no file and silently rendered nothing.
    if uploaded_file is None and st.session_state.sample_loaded:
        uploaded_file = SampleFile(SAMPLE_FASTA, "sample_sequences.fasta")

    # An explicit upload takes over from the sample.
    if uploaded_file is not None and not isinstance(uploaded_file, SampleFile):
        st.session_state.sample_loaded = False

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
    st.markdown(
        section_header(
            'ruler', 'Analysis Results',
            'Computed directly from the uploaded sequence(s) — no external calls or model in this path.'
        ),
        unsafe_allow_html=True,
    )

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

    # Two per row rather than three: at sidebar + narrow-viewport widths a
    # three-way split truncates every chart title, and this row already
    # matches the width the composition/length pair was designed at.
    sequences = st.session_state.sequences or []
    col1, col2 = st.columns(2)

    with col1:
        composition = results.get('nucleotide_composition', {})
        if composition:
            fig = create_composition_plot(composition)
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        analyzer = st.session_state.analyzer
        gc_contents = [analyzer.calculate_gc_content(str(seq.seq)) for seq in sequences]
        labels = [seq.id for seq in sequences]
        if gc_contents:
            fig = create_gc_plot(gc_contents, labels)
            # A per-sequence legend is useful at a handful of records; past
            # that it turns into unreadable colour soup, so it's dropped
            # rather than left to overrun the plot.
            fig.update_layout(showlegend=len(labels) <= 12, legend_title_text="Sequence")
            st.plotly_chart(fig, use_container_width=True)

    lengths = [len(seq.seq) for seq in sequences]
    if lengths:
        fig = create_length_distribution(lengths)
        st.plotly_chart(fig, use_container_width=True)

    # Quality scores only exist for FASTQ input — Biopython attaches
    # per-base Phred scores as letter_annotations only when parsing that
    # format, so this section has nothing to show for FASTA and stays hidden
    # rather than rendering an empty chart.
    if results.get('format') == 'fastq':
        analyzer = st.session_state.analyzer
        quality_stats = analyzer.calculate_quality_statistics(sequences)
        if quality_stats:
            st.markdown(
                section_header(
                    'ruler', 'Quality Scores',
                    'Per-base Phred quality from the FASTQ file, pooled across all reads.'
                ),
                unsafe_allow_html=True,
            )
            q1, q2, q3, q4 = st.columns(4)
            q1.metric("Mean Q", f"{quality_stats['average_quality']:.1f}")
            q2.metric("Median Q", f"{quality_stats['median_quality']:.1f}")
            total = quality_stats['total_bases']
            q3.metric("Bases ≥ Q20", f"{quality_stats['q20_bases'] / total * 100:.1f}%" if total else "—")
            q4.metric("Bases ≥ Q30", f"{quality_stats['q30_bases'] / total * 100:.1f}%" if total else "—")

            quality_data = [
                seq.letter_annotations['phred_quality']
                for seq in sequences
                if 'phred_quality' in seq.letter_annotations
            ]
            if quality_data:
                fig = create_quality_plot(quality_data)
                st.plotly_chart(fig, use_container_width=True)

    # Feature annotations only exist for GenBank input — Biopython attaches
    # `.features` only when parsing that format, so FASTA/FASTQ have nothing
    # to show here and this section stays hidden for them.
    if results.get('format') == 'genbank':
        render_genbank_features(sequences)


def render_genbank_features(sequences: list) -> None:
    """Render the feature annotations carried by parsed GenBank records."""
    rows = []
    for seq in sequences[:10]:
        for feature in extract_genbank_features(seq):
            qualifiers = feature['qualifiers']
            label = (qualifiers.get('gene') or qualifiers.get('product') or [''])[0]
            rows.append({
                'Sequence': seq.id,
                'Type': feature['type'],
                'Start': feature['start'],
                'End': feature['end'],
                'Strand': feature['strand'] or '',
                'Gene / Product': label,
            })

    st.markdown(
        section_header(
            'sequenceFile', 'GenBank Feature Annotations',
            'Features (source, CDS, gene, etc.) as recorded in the uploaded GenBank file — not computed.'
        ),
        unsafe_allow_html=True,
    )

    if not rows:
        st.info("This GenBank file has no FEATURES entries to show.")
        return

    import pandas as pd
    st.dataframe(pd.DataFrame(rows[:200]), use_container_width=True, hide_index=True)
    if len(rows) > 200:
        st.caption(f"Showing first 200 of {len(rows)} features.")


def render_sidebar():
    """Render the sidebar with analysis options."""
    with st.sidebar:
        st.markdown("### Analysis Options")
        
        st.markdown("---")
        
        sidebar_section_header('kmer', 'K-mer Analysis')
        kmer_size = st.slider("K-mer size", min_value=2, max_value=10, value=3)
        run_kmer = st.button("Run K-mer Analysis")

        st.markdown("---")

        sidebar_section_header('readingFrame', 'ORF Detection')
        min_orf_length = st.number_input("Minimum ORF length (bp)",
                                          min_value=30, max_value=1000, value=100)
        include_reverse = st.checkbox("Include reverse complement", value=True)
        run_orf = st.button("Find ORFs")

        st.markdown("---")

        # Reuses the ORF Detection settings above rather than duplicating a
        # second min-length/reverse-complement control: the per-ORF codon
        # table is built from exactly the ORFs those settings would find.
        sidebar_section_header('readingFrame', 'Codon Usage')
        st.caption("Uses the ORF Detection settings above to find ORFs.")
        run_codon = st.button("Analyze Codon Usage")

        st.markdown("---")

        sidebar_section_header('motif', 'Motif Search')
        motif_pattern = st.text_input("Pattern (IUPAC)", placeholder="ATGR")
        enzyme = st.selectbox("Or select enzyme",
                              ["Custom"] + list(get_enzyme_list().keys()))
        run_motif = st.button("Search Motif")

        st.markdown("---")

        # Alignment is the one analysis that needs TWO sequences, so it gets its
        # own pair of selectors rather than running over the loaded set.
        sidebar_section_header('alignment', 'Pairwise Alignment')
        records = st.session_state.sequences or []
        labels = [f"{i + 1}. {r.id}" for i, r in enumerate(records)]

        if len(records) < 2:
            st.caption(
                "Needs two sequences. Upload a FASTA containing at least two "
                "records to enable alignment."
            )
            align_i = align_j = None
            align_mode = "global"
            run_align = False
        else:
            align_i = st.selectbox("Sequence A", labels, index=0, key="align_a")
            align_j = st.selectbox("Sequence B", labels, index=1, key="align_b")
            align_mode = st.radio(
                "Mode",
                ["global", "local"],
                horizontal=True,
                help=(
                    "Global (Needleman-Wunsch) aligns the sequences end to end. "
                    "Local (Smith-Waterman) finds the best-scoring subregion."
                ),
            )
            run_align = st.button("Align Sequences")

        st.markdown("---")

        sidebar_section_header('sequenceFile', 'Export')
        export_format = st.selectbox("Format", ["CSV", "JSON", "ZIP"])
        run_export = st.button("Export Results")

        return {
            'kmer': {'run': run_kmer, 'size': kmer_size},
            'orf': {'run': run_orf, 'min_length': min_orf_length,
                    'reverse': include_reverse},
            'codon': {'run': run_codon, 'min_length': min_orf_length,
                      'reverse': include_reverse},
            'motif': {'run': run_motif, 'pattern': motif_pattern,
                      'enzyme': enzyme},
            'align': {'run': run_align, 'a': align_i, 'b': align_j,
                      'mode': align_mode, 'labels': labels},
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
        
        st.markdown(
            section_header('kmer', f'{size}-mer Analysis',
                           f'Frequency of every {size}-base window across the loaded sequence(s), top 20 shown.'),
            unsafe_allow_html=True,
        )

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
        
        reverse_note = " and its reverse complement" if include_reverse else ""
        st.markdown(
            section_header('readingFrame', 'Open Reading Frames',
                           f'Start-to-stop spans of at least {min_length} bp across all 3 forward frames{reverse_note}.'),
            unsafe_allow_html=True,
        )

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


def _render_codon_table(usage) -> None:
    """Render one CodonUsageResult as a table, showing only codons whose
    amino-acid family was actually observed. Rows for codons that were never
    used within an observed family show a real 0 count and RSCU=0.0 (that's
    data); families with zero observations are omitted entirely rather than
    padded with a fabricated 0 (that would look like data but isn't)."""
    import pandas as pd

    rows = [r for r in get_codon_table_rows(usage) if r['fraction_of_aa_pct'] is not None]

    if not rows:
        st.info("No standard-base codons were counted for this table.")
    else:
        df = pd.DataFrame(rows).rename(columns={
            'codon': 'Codon', 'amino_acid': 'AA', 'count': 'Count',
            'fraction_of_aa_pct': '% of AA', 'rscu': 'RSCU',
        })
        df['% of AA'] = df['% of AA'].round(1)
        df['RSCU'] = df['RSCU'].round(2)
        st.dataframe(df, use_container_width=True, hide_index=True)

    # Ambiguous and trailing-partial codons are never folded into the table
    # above or silently dropped — always state the count, even when it's 0.
    st.caption(
        f"{usage.counted_codons:,} standard codons counted &middot; "
        f"{usage.ambiguous_codons:,} ambiguous (non-ACGT) &middot; "
        f"{usage.partial_trailing_bases:,} trailing base(s) outside a full codon."
    )


def run_codon_usage_analysis(min_length: int, include_reverse: bool):
    """Run codon usage analysis and display per-ORF and whole-sequence tables."""
    sequences = st.session_state.sequences
    if not sequences:
        return

    with st.spinner("Counting codons..."):
        analyzed = sequences[:10]

        # Whole sequence: a raw frame-1 tally of the sequence as uploaded,
        # not filtered to coding regions.
        whole_seq_usage = merge_codon_usage(
            [count_codons(str(seq.seq), frame=0) for seq in analyzed]
        )

        # Per-ORF: every ORF the current ORF Detection settings find across
        # the same sequences, each counted from its own start codon, then
        # aggregated into one table.
        all_orfs = []
        for seq in analyzed:
            all_orfs.extend(find_orfs(
                str(seq.seq), min_length=min_length, include_reverse=include_reverse
            ))
        orf_usage = merge_codon_usage([codon_usage_for_orf(orf) for orf in all_orfs]) if all_orfs else None

        st.markdown(
            section_header(
                'readingFrame', 'Codon Usage',
                'Codon counts and relative synonymous codon usage (RSCU) — '
                'RSCU=1.0 is even use across an amino acid\'s synonymous codons, '
                '>1.0 favored, <1.0 avoided.'
            ),
            unsafe_allow_html=True,
        )

        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f"**Whole sequence** &middot; {len(analyzed)} sequence(s), frame 1")
            _render_codon_table(whole_seq_usage)
        with col2:
            if orf_usage is not None:
                st.markdown(f"**Detected ORFs** &middot; aggregate over {len(all_orfs)} ORF(s)")
                _render_codon_table(orf_usage)
            else:
                st.markdown("**Detected ORFs**")
                st.info(
                    f"No ORFs found at the current ORF Detection minimum "
                    f"length ({min_length} bp), so there is no per-ORF codon "
                    f"table. The whole-sequence tally is unaffected."
                )


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
        
        st.markdown(
            section_header('motif', f'Motif Search: {pattern}',
                           'IUPAC-ambiguous pattern match across the loaded sequence(s).'),
            unsafe_allow_html=True,
        )
        display_motif_results(matches)


def run_alignment(label_a: str, label_b: str, labels: list, mode: str):
    """Align two of the loaded sequences and display the result."""
    records = st.session_state.sequences or []
    if len(records) < 2:
        st.warning("Alignment needs at least two sequences.")
        return

    try:
        i, j = labels.index(label_a), labels.index(label_b)
    except ValueError:
        st.warning("Select two sequences to align.")
        return

    if i == j:
        st.warning("Pick two different sequences — aligning a sequence to itself is trivially 100% identical.")
        return

    rec_a, rec_b = records[i], records[j]
    with st.spinner("Aligning..."):
        result = align_sequences(str(rec_a.seq), str(rec_b.seq), mode=mode)

    st.markdown(section_header('alignment', 'Pairwise Alignment'), unsafe_allow_html=True)
    st.caption(
        f"{rec_a.id} ({len(rec_a.seq):,} bp) vs {rec_b.id} ({len(rec_b.seq):,} bp) · "
        f"{mode} alignment (Biopython PairwiseAligner)"
    )

    # A failed alignment used to render as "Identity 0.0% · Score 0.0 · Gaps 0",
    # which is exactly what a real comparison of two unrelated sequences looks
    # like. Report the failure instead of four numbers that describe nothing.
    if not result.ok:
        st.error(
            f"**The alignment did not run, so there are no results to show.**\n\n"
            f"Reason: {result.error}\n\n"
            f"This is a failure, not a finding — it does not mean the two "
            f"sequences are dissimilar."
        )
        return

    # Distinct from the failure above: the aligner ran and found nothing. That
    # is a real answer about these two sequences, so say so rather than
    # printing 0.0% across four metrics as though it were a measurement.
    if result.alignment_length == 0:
        st.info(
            f"**No {mode} alignment exists between these two sequences.**\n\n"
            f"The aligner ran successfully and found no alignable region — "
            f"this is a result, not an error."
        )
        return

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Identity", f"{result.identity:.1f}%",
              help=f"{result.identity_count:,} of {result.alignment_length:,} aligned positions match")
    c2.metric("Score", f"{result.score:,.1f}")
    c3.metric("Gaps", f"{result.gaps:,}")
    c4.metric("Aligned length", f"{result.alignment_length:,}")

    # Monospace, fixed width: the match line only lines up in a mono font.
    st.code(format_alignment_display(result), language=None)

    st.download_button(
        "Download alignment (.txt)",
        data=format_alignment_display(result),
        file_name=f"alignment_{rec_a.id}_vs_{rec_b.id}.txt".replace("/", "_"),
        mime="text/plain",
    )


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
            GenomeSight v1.0.0
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

            if sidebar_options['codon']['run']:
                run_codon_usage_analysis(
                    sidebar_options['codon']['min_length'],
                    sidebar_options['codon']['reverse']
                )

            if sidebar_options['motif']['run']:
                run_motif_analysis(
                    sidebar_options['motif']['pattern'],
                    sidebar_options['motif']['enzyme']
                )
            
            if sidebar_options['align']['run']:
                run_alignment(
                    sidebar_options['align']['a'],
                    sidebar_options['align']['b'],
                    sidebar_options['align']['labels'],
                    sidebar_options['align']['mode'],
                )

            if sidebar_options['export']['run']:
                run_export(sidebar_options['export']['format'])
    
    render_footer()


if __name__ == "__main__":
    main()
