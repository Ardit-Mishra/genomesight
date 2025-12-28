"""
Core Sequence Analysis Module

This module provides the main SequenceAnalyzer class for DNA/RNA sequence analysis.
It includes functions for GC content calculation, nucleotide composition,
complexity analysis, and comprehensive sequence statistics.

Design Notes:
    - All functions are pure Python with no Streamlit dependencies
    - Uses Biopython for sequence handling
    - NumPy for numerical computations
    - Type hints throughout for code clarity

Main Classes:
    SequenceAnalyzer: Primary analysis class with all sequence metrics
"""

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from collections import Counter
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class SequenceStats:
    """Data class for sequence statistics."""
    sequence_count: int
    total_length: int
    average_length: float
    median_length: float
    min_length: int
    max_length: int
    average_gc: float
    composition: Dict[str, int]


class SequenceAnalyzer:
    """
    Comprehensive sequence analyzer for genomic data.
    
    This class provides methods for analyzing DNA and RNA sequences,
    including GC content, nucleotide composition, complexity metrics,
    and k-mer analysis.
    
    Attributes:
        valid_nucleotides: Set of valid DNA nucleotide characters
        complement_map: Dictionary mapping nucleotides to their complements
    
    Example:
        >>> analyzer = SequenceAnalyzer()
        >>> gc = analyzer.calculate_gc_content("ATGCATGC")
        >>> print(f"GC Content: {gc:.2f}%")
    """
    
    def __init__(self):
        """Initialize the SequenceAnalyzer with nucleotide mappings."""
        self.valid_nucleotides = set('ATCGN')
        self.complement_map = {
            'A': 'T', 'T': 'A', 
            'C': 'G', 'G': 'C', 
            'N': 'N'
        }
        self.codon_table = self._build_codon_table()
    
    def _build_codon_table(self) -> Dict[str, str]:
        """Build the standard genetic code translation table."""
        return {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
    
    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate the GC content percentage of a DNA sequence.
        
        GC content is the proportion of guanine (G) and cytosine (C) bases
        in a nucleotide sequence. This metric is important for:
        - Species identification
        - Thermal stability prediction
        - PCR primer design
        
        Args:
            sequence: DNA sequence string (case-insensitive)
        
        Returns:
            GC content as a percentage (0-100)
        
        Example:
            >>> analyzer.calculate_gc_content("ATGCATGC")
            50.0
        """
        if not sequence:
            return 0.0
        
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        total_count = len([n for n in sequence if n in 'ATCG'])
        
        if total_count == 0:
            return 0.0
        
        return (gc_count / total_count) * 100
    
    def calculate_nucleotide_composition(self, sequence: str) -> Dict[str, int]:
        """
        Count the occurrence of each nucleotide in a sequence.
        
        Args:
            sequence: DNA/RNA sequence string
        
        Returns:
            Dictionary with nucleotide counts {'A': n, 'T': n, 'C': n, 'G': n, 'N': n}
        """
        sequence = sequence.upper()
        composition = Counter(sequence)
        
        for nucleotide in 'ATCGN':
            if nucleotide not in composition:
                composition[nucleotide] = 0
        
        return dict(composition)
    
    def reverse_complement(self, sequence: str) -> str:
        """
        Generate the reverse complement of a DNA sequence.
        
        The reverse complement is found by complementing each base
        and then reversing the entire sequence.
        
        Args:
            sequence: DNA sequence string
        
        Returns:
            Reverse complement sequence
        
        Example:
            >>> analyzer.reverse_complement("ATGC")
            "GCAT"
        """
        sequence = sequence.upper()
        complement = ''.join(
            self.complement_map.get(base, base) for base in sequence
        )
        return complement[::-1]
    
    def calculate_complexity(self, sequence: str, window_size: int = 10) -> float:
        """
        Calculate sequence complexity using Shannon entropy.
        
        Higher entropy indicates more complex/random sequences.
        Lower entropy indicates repetitive or low-complexity regions.
        
        Args:
            sequence: DNA sequence string
            window_size: Size of sliding window for calculation
        
        Returns:
            Average Shannon entropy (bits)
        """
        if len(sequence) < window_size:
            return 0.0
        
        complexities = []
        
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size].upper()
            composition = Counter(window)
            total = len(window)
            
            entropy = 0.0
            for count in composition.values():
                if count > 0:
                    p = count / total
                    entropy -= p * np.log2(p)
            
            complexities.append(entropy)
        
        return np.mean(complexities) if complexities else 0.0
    
    def calculate_sequence_statistics(self, sequences: List) -> Dict[str, Any]:
        """
        Calculate comprehensive statistics for a list of sequences.
        
        Args:
            sequences: List of Bio.SeqRecord objects
        
        Returns:
            Dictionary containing all computed statistics
        """
        if not sequences:
            return {}
        
        lengths = [len(seq.seq) for seq in sequences]
        gc_contents = [self.calculate_gc_content(str(seq.seq)) for seq in sequences]
        
        total_composition = Counter()
        for seq in sequences:
            seq_comp = self.calculate_nucleotide_composition(str(seq.seq))
            for nucleotide, count in seq_comp.items():
                total_composition[nucleotide] += count
        
        total_bases = sum(total_composition.values())
        composition_pct = {
            nuc: (count / total_bases) * 100 
            for nuc, count in total_composition.items()
        } if total_bases > 0 else {}
        
        return {
            'sequence_count': len(sequences),
            'total_length': sum(lengths),
            'average_length': float(np.mean(lengths)),
            'median_length': float(np.median(lengths)),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'length_std': float(np.std(lengths)),
            'average_gc_content': float(np.mean(gc_contents)),
            'median_gc_content': float(np.median(gc_contents)),
            'min_gc_content': float(min(gc_contents)),
            'max_gc_content': float(max(gc_contents)),
            'nucleotide_composition': dict(total_composition),
            'composition_percentages': composition_pct
        }
    
    def calculate_quality_statistics(self, sequences: List) -> Optional[Dict[str, float]]:
        """
        Calculate quality score statistics for FASTQ sequences.
        
        Args:
            sequences: List of Bio.SeqRecord objects with quality annotations
        
        Returns:
            Dictionary of quality metrics, or None if no quality data
        """
        try:
            quality_scores = []
            for seq in sequences:
                if hasattr(seq, 'letter_annotations'):
                    if 'phred_quality' in seq.letter_annotations:
                        quality_scores.extend(seq.letter_annotations['phred_quality'])
            
            if not quality_scores:
                return None
            
            return {
                'average_quality': float(np.mean(quality_scores)),
                'median_quality': float(np.median(quality_scores)),
                'min_quality': min(quality_scores),
                'max_quality': max(quality_scores),
                'quality_std': float(np.std(quality_scores)),
                'q20_bases': sum(1 for q in quality_scores if q >= 20),
                'q30_bases': sum(1 for q in quality_scores if q >= 30),
                'total_bases': len(quality_scores)
            }
        except Exception as e:
            logger.warning(f"Quality calculation failed: {e}")
            return None
    
    def analyze_kmers(self, sequences: List, k: int = 3) -> Dict[str, Any]:
        """
        Analyze k-mer frequencies in sequences.
        
        K-mers are subsequences of length k. This analysis is useful for:
        - Codon usage analysis (k=3)
        - Motif discovery
        - Sequence comparison
        - Repeat identification
        
        Args:
            sequences: List of Bio.SeqRecord objects
            k: K-mer length (default 3 for codons)
        
        Returns:
            Dictionary with k-mer statistics and counts
        """
        all_kmers = Counter()
        
        for seq_record in sequences:
            sequence = str(seq_record.seq).upper()
            
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                if all(base in 'ATCG' for base in kmer):
                    all_kmers[kmer] += 1
        
        total = sum(all_kmers.values())
        
        return {
            'k': k,
            'total_kmers': total,
            'unique_kmers': len(all_kmers),
            'most_common': all_kmers.most_common(20),
            'all_kmers': dict(all_kmers)
        }
    
    def translate_codon(self, codon: str) -> str:
        """
        Translate a single codon to its amino acid.
        
        Args:
            codon: Three-letter DNA codon
        
        Returns:
            Single-letter amino acid code, or 'X' for unknown
        """
        return self.codon_table.get(codon.upper(), 'X')
