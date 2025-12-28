"""
Export and Report Generation Module

This module provides functions for exporting analysis results in various
formats including JSON, CSV, and comprehensive text reports. It also
handles ZIP archive creation for batch downloads.

Export Formats:
    - JSON: Structured data for API consumption
    - CSV: Spreadsheet-compatible tabular data
    - Text: Human-readable analysis reports
    - ZIP: Archive containing all export formats
"""

import json
import zipfile
from io import BytesIO
from datetime import datetime
from typing import Dict, Any, List, Optional
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def generate_json_report(results: Dict[str, Any]) -> str:
    """
    Generate JSON format analysis report.
    
    Args:
        results: Analysis results dictionary
    
    Returns:
        JSON formatted string
    """
    report = {
        'metadata': {
            'generated_at': datetime.now().isoformat(),
            'tool': 'Genome Sequencing Analyzer',
            'version': '1.0.0'
        },
        'results': _serialize_results(results)
    }
    
    return json.dumps(report, indent=2, default=str)


def _serialize_results(obj: Any) -> Any:
    """Recursively serialize objects for JSON."""
    if hasattr(obj, '__dict__'):
        return {k: _serialize_results(v) for k, v in obj.__dict__.items()}
    elif isinstance(obj, dict):
        return {k: _serialize_results(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_serialize_results(item) for item in obj]
    elif hasattr(obj, 'tolist'):
        return obj.tolist()
    else:
        return obj


def generate_csv_report(results: Dict[str, Dict]) -> str:
    """
    Generate CSV format analysis report.
    
    Args:
        results: Analysis results by filename
    
    Returns:
        CSV formatted string
    """
    rows = []
    
    for filename, data in results.items():
        row = {
            'Filename': filename,
            'Format': data.get('format', 'unknown').upper(),
            'Sequence_Count': data.get('sequence_count', 0),
            'Total_Length': data.get('total_length', 0),
            'Average_Length': round(data.get('average_length', 0), 2),
            'Min_Length': data.get('min_length', 0),
            'Max_Length': data.get('max_length', 0),
            'Average_GC': round(data.get('average_gc_content', 0), 2)
        }
        
        composition = data.get('nucleotide_composition', {})
        for nuc in 'ATCGN':
            row[f'{nuc}_Count'] = composition.get(nuc, 0)
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    return df.to_csv(index=False)


def generate_text_report(results: Dict[str, Dict]) -> str:
    """
    Generate comprehensive text analysis report.
    
    Args:
        results: Analysis results by filename
    
    Returns:
        Formatted text report string
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    lines = [
        "=" * 70,
        "GENOME SEQUENCING ANALYSIS REPORT",
        "=" * 70,
        f"Generated: {timestamp}",
        f"Tool: Genome Sequencing Analyzer v1.0.0",
        "",
        "SUMMARY",
        "-" * 40,
        f"Total Files Analyzed: {len(results)}",
        f"Total Sequences: {sum(d.get('sequence_count', 0) for d in results.values()):,}",
        f"Total Nucleotides: {sum(d.get('total_length', 0) for d in results.values()):,}",
        "",
        "DETAILED RESULTS",
        "-" * 40
    ]
    
    for filename, data in results.items():
        lines.extend([
            "",
            f"File: {filename}",
            f"  Format: {data.get('format', 'unknown').upper()}",
            f"  Sequences: {data.get('sequence_count', 0):,}",
            f"  Total Length: {data.get('total_length', 0):,} bp",
            f"  Average Length: {data.get('average_length', 0):.2f} bp",
            f"  Length Range: {data.get('min_length', 0)} - {data.get('max_length', 0)} bp",
            f"  Average GC Content: {data.get('average_gc_content', 0):.2f}%"
        ])
        
        composition = data.get('nucleotide_composition', {})
        total = data.get('total_length', 1)
        lines.append("")
        lines.append("  Nucleotide Composition:")
        for nuc in 'ATCG':
            count = composition.get(nuc, 0)
            pct = (count / total * 100) if total > 0 else 0
            lines.append(f"    {nuc}: {count:,} ({pct:.1f}%)")
    
    lines.extend([
        "",
        "=" * 70,
        "End of Report",
        "=" * 70
    ])
    
    return '\n'.join(lines)


def create_download_zip(results: Dict[str, Dict]) -> bytes:
    """
    Create ZIP archive containing all export formats.
    
    Args:
        results: Analysis results dictionary
    
    Returns:
        ZIP file content as bytes
    """
    zip_buffer = BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        zf.writestr('analysis_report.csv', generate_csv_report(results))
        zf.writestr('analysis_report.json', generate_json_report(results))
        zf.writestr('analysis_report.txt', generate_text_report(results))
        
        readme = """GENOME ANALYSIS RESULTS
=======================

This archive contains your genome sequence analysis results:

- analysis_report.csv  : Spreadsheet-ready data
- analysis_report.json : Structured data for APIs  
- analysis_report.txt  : Human-readable report

Generated by Genome Sequencing Analyzer v1.0.0
"""
        zf.writestr('README.txt', readme)
    
    zip_buffer.seek(0)
    logger.info("Created download ZIP archive")
    return zip_buffer.getvalue()


def export_orfs_csv(orfs: List) -> str:
    """
    Export ORF data to CSV format.
    
    Args:
        orfs: List of ORF objects
    
    Returns:
        CSV formatted string
    """
    rows = []
    for i, orf in enumerate(orfs, 1):
        rows.append({
            'ORF_ID': i,
            'Start': orf.start,
            'End': orf.end,
            'Length_nt': orf.length_nt,
            'Length_aa': orf.length_aa,
            'Frame': orf.frame,
            'Strand': orf.strand,
            'Start_Codon': orf.start_codon,
            'Stop_Codon': orf.stop_codon,
            'GC_Content': round(orf.gc_content, 2),
            'Protein': orf.protein[:50] + '...' if len(orf.protein) > 50 else orf.protein
        })
    
    df = pd.DataFrame(rows)
    return df.to_csv(index=False)


def export_motif_matches_csv(matches: List) -> str:
    """
    Export motif match data to CSV format.
    
    Args:
        matches: List of MotifMatch objects
    
    Returns:
        CSV formatted string
    """
    rows = []
    for match in matches:
        rows.append({
            'Pattern': match.pattern,
            'Sequence_ID': match.sequence_id,
            'Start_0based': match.start_0based,
            'Start_1based': match.start_1based,
            'End': match.end,
            'Matched': match.matched_sequence,
            'Context': match.context
        })
    
    df = pd.DataFrame(rows)
    return df.to_csv(index=False)
