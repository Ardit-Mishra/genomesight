import streamlit as st
from Bio import SeqIO
from io import StringIO
import pandas as pd
import json
from datetime import datetime
from typing import Dict, List, Any, Optional, Union

def validate_file_format(file_content: str, filename: str) -> Optional[str]:
    """
    Validate file format based on content and filename
    
    Args:
        file_content (str): Content of the uploaded file
        filename (str): Name of the uploaded file
        
    Returns:
        Optional[str]: File format ('fasta' or 'fastq') or None if invalid
    """
    try:
        # Check file extension
        if filename.lower().endswith(('.fasta', '.fa')):
            # Validate FASTA format
            if file_content.startswith('>'):
                return 'fasta'
        elif filename.lower().endswith(('.fastq', '.fq')):
            # Validate FASTQ format
            lines = file_content.strip().split('\n')
            if len(lines) >= 4 and lines[0].startswith('@') and lines[2].startswith('+'):
                return 'fastq'
        
        # Try to auto-detect format based on content
        if file_content.startswith('>'):
            return 'fasta'
        elif file_content.startswith('@'):
            return 'fastq'
        
        return None
    except Exception:
        return None

def parse_uploaded_files(uploaded_files) -> Dict[str, Dict[str, Any]]:
    """
    Parse uploaded sequence files
    
    Args:
        uploaded_files: Streamlit uploaded files
        
    Returns:
        Dict[str, Dict[str, Any]]: Parsed file data
    """
    parsed_files = {}
    
    for uploaded_file in uploaded_files:
        try:
            # Read file content
            content = uploaded_file.read()
            
            # Try to decode as text
            try:
                if isinstance(content, bytes):
                    content = content.decode('utf-8')
            except UnicodeDecodeError:
                st.error(f"Could not decode file {uploaded_file.name}. Please ensure it's a text file.")
                continue
            
            # Validate format
            file_format = validate_file_format(content, uploaded_file.name)
            
            if not file_format:
                st.warning(f"Invalid format for file {uploaded_file.name}. Skipping.")
                continue
            
            # Parse sequences
            content_io = StringIO(content)
            sequences = list(SeqIO.parse(content_io, file_format))
            
            if not sequences:
                st.warning(f"No sequences found in {uploaded_file.name}. Skipping.")
                continue
            
            parsed_files[uploaded_file.name] = {
                'format': file_format,
                'sequences': sequences,
                'content': content
            }
            
        except Exception as e:
            st.error(f"Error parsing file {uploaded_file.name}: {str(e)}")
            continue
    
    return parsed_files

def generate_report(results: Dict[str, Dict[str, Any]], format_type: str) -> Union[str, bytes]:
    """
    Generate analysis report in specified format
    
    Args:
        results (Dict[str, Dict[str, Any]]): Analysis results
        format_type (str): Output format ('csv', 'json', 'full')
        
    Returns:
        Union[str, bytes]: Generated report content
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    if format_type == 'csv':
        # Create CSV report
        rows = []
        for filename, data in results.items():
            row = {
                'Filename': filename,
                'Format': data['format'].upper(),
                'Sequence_Count': data['sequence_count'],
                'Total_Length': data['total_length'],
                'Average_Length': round(data['average_length'], 2),
                'Min_Length': data['min_length'],
                'Max_Length': data['max_length'],
                'Average_GC_Content': round(data['average_gc_content'], 2),
                'A_Count': data['nucleotide_composition'].get('A', 0),
                'T_Count': data['nucleotide_composition'].get('T', 0),
                'C_Count': data['nucleotide_composition'].get('C', 0),
                'G_Count': data['nucleotide_composition'].get('G', 0),
                'N_Count': data['nucleotide_composition'].get('N', 0)
            }
            
            # Add quality statistics if available
            if data.get('quality_statistics'):
                quality_stats = data['quality_statistics']
                row.update({
                    'Average_Quality': round(quality_stats['average_quality'], 2),
                    'Min_Quality': quality_stats['min_quality'],
                    'Max_Quality': quality_stats['max_quality'],
                    'Q20_Bases': quality_stats['q20_bases'],
                    'Q30_Bases': quality_stats['q30_bases']
                })
            
            rows.append(row)
        
        df = pd.DataFrame(rows)
        return df.to_csv(index=False)
    
    elif format_type == 'json':
        # Create JSON report
        report = {
            'metadata': {
                'generated_at': timestamp,
                'tool': 'Genome Sequencing Analyzer',
                'version': '1.0.0'
            },
            'summary': {
                'total_files': len(results),
                'total_sequences': sum(data['sequence_count'] for data in results.values()),
                'total_nucleotides': sum(data['total_length'] for data in results.values())
            },
            'results': results
        }
        
        return json.dumps(report, indent=2, default=str)
    
    elif format_type == 'full':
        # Create comprehensive text report
        report_lines = [
            "=" * 80,
            "GENOME SEQUENCING ANALYSIS REPORT",
            "=" * 80,
            f"Generated: {timestamp}",
            f"Tool: Genome Sequencing Analyzer v1.0.0",
            "",
            "SUMMARY",
            "-" * 40,
            f"Total Files Analyzed: {len(results)}",
            f"Total Sequences: {sum(data['sequence_count'] for data in results.values()):,}",
            f"Total Nucleotides: {sum(data['total_length'] for data in results.values()):,}",
            "",
            "DETAILED RESULTS",
            "-" * 40
        ]
        
        for filename, data in results.items():
            report_lines.extend([
                "",
                f"File: {filename}",
                f"  Format: {data['format'].upper()}",
                f"  Sequences: {data['sequence_count']:,}",
                f"  Total Length: {data['total_length']:,} bp",
                f"  Average Length: {data['average_length']:.2f} bp",
                f"  Length Range: {data['min_length']} - {data['max_length']} bp",
                f"  Average GC Content: {data['average_gc_content']:.2f}%",
                "",
                "  Nucleotide Composition:",
                f"    A: {data['nucleotide_composition'].get('A', 0):,} ({data['nucleotide_composition'].get('A', 0)/data['total_length']*100:.1f}%)",
                f"    T: {data['nucleotide_composition'].get('T', 0):,} ({data['nucleotide_composition'].get('T', 0)/data['total_length']*100:.1f}%)",
                f"    C: {data['nucleotide_composition'].get('C', 0):,} ({data['nucleotide_composition'].get('C', 0)/data['total_length']*100:.1f}%)",
                f"    G: {data['nucleotide_composition'].get('G', 0):,} ({data['nucleotide_composition'].get('G', 0)/data['total_length']*100:.1f}%)",
                f"    N: {data['nucleotide_composition'].get('N', 0):,} ({data['nucleotide_composition'].get('N', 0)/data['total_length']*100:.1f}%)"
            ])
            
            # Add quality statistics if available
            if data.get('quality_statistics'):
                quality_stats = data['quality_statistics']
                report_lines.extend([
                    "",
                    "  Quality Statistics:",
                    f"    Average Quality: {quality_stats['average_quality']:.2f}",
                    f"    Quality Range: {quality_stats['min_quality']} - {quality_stats['max_quality']}",
                    f"    Q20 Bases: {quality_stats['q20_bases']:,} ({quality_stats['q20_bases']/quality_stats['total_bases']*100:.1f}%)",
                    f"    Q30 Bases: {quality_stats['q30_bases']:,} ({quality_stats['q30_bases']/quality_stats['total_bases']*100:.1f}%)"
                ])
        
        report_lines.extend([
            "",
            "=" * 80,
            "End of Report",
            "=" * 80
        ])
        
        return '\n'.join(report_lines)
    
    return ""

def format_large_number(number: int) -> str:
    """
    Format large numbers with appropriate units
    
    Args:
        number (int): Number to format
        
    Returns:
        str: Formatted number string
    """
    if number >= 1_000_000_000:
        return f"{number / 1_000_000_000:.1f}B"
    elif number >= 1_000_000:
        return f"{number / 1_000_000:.1f}M"
    elif number >= 1_000:
        return f"{number / 1_000:.1f}K"
    else:
        return str(number)

def calculate_file_checksum(content: str) -> str:
    """
    Calculate MD5 checksum for file content
    
    Args:
        content (str): File content
        
    Returns:
        str: MD5 checksum
    """
    import hashlib
    return hashlib.md5(content.encode()).hexdigest()

def validate_sequence_characters(sequence: str) -> bool:
    """
    Validate that sequence contains only valid nucleotide characters
    
    Args:
        sequence (str): DNA/RNA sequence
        
    Returns:
        bool: True if valid, False otherwise
    """
    valid_chars = set('ATCGRYSWKMBDHVN')  # Standard IUPAC nucleotide codes
    return all(char.upper() in valid_chars for char in sequence)

def get_sequence_type(sequence: str) -> str:
    """
    Determine if sequence is DNA or RNA based on nucleotide composition
    
    Args:
        sequence (str): Nucleotide sequence
        
    Returns:
        str: 'DNA', 'RNA', or 'Unknown'
    """
    sequence = sequence.upper()
    has_t = 'T' in sequence
    has_u = 'U' in sequence
    
    if has_u and not has_t:
        return 'RNA'
    elif has_t and not has_u:
        return 'DNA'
    elif has_t and has_u:
        return 'Mixed (DNA/RNA)'
    else:
        return 'Unknown'

@st.cache_data
def load_sample_data():
    """
    Load sample genomic data for demonstration purposes
    Note: This function is cached for performance
    
    Returns:
        Dict: Sample data structure
    """
    # This would normally load from a file or database
    # For now, return empty dict since we don't want mock data
    return {}

def export_sequences_to_fasta(sequences: List, filename: str) -> str:
    """
    Export sequences to FASTA format
    
    Args:
        sequences (List): List of Bio.SeqRecord objects
        filename (str): Output filename
        
    Returns:
        str: FASTA formatted string
    """
    fasta_lines = []
    
    for i, seq in enumerate(sequences):
        seq_id = seq.id if hasattr(seq, 'id') else f"sequence_{i+1}"
        description = seq.description if hasattr(seq, 'description') else ""
        
        # Create header line
        header = f">{seq_id}"
        if description and description != seq_id:
            header += f" {description}"
        
        fasta_lines.append(header)
        
        # Add sequence in 80-character lines
        sequence_str = str(seq.seq)
        for j in range(0, len(sequence_str), 80):
            fasta_lines.append(sequence_str[j:j+80])
    
    return '\n'.join(fasta_lines)
