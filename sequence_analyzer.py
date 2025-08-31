from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from collections import Counter
import re
from typing import List, Dict, Any, Optional, Tuple

class SequenceAnalyzer:
    """
    A comprehensive sequence analyzer for genomic data analysis
    """
    
    def __init__(self):
        self.valid_nucleotides = set('ATCGN')
        self.complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    
    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate GC content of a sequence
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            float: GC content percentage
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
        Calculate nucleotide composition of a sequence
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            Dict[str, int]: Dictionary with nucleotide counts
        """
        sequence = sequence.upper()
        composition = Counter(sequence)
        
        # Ensure all standard nucleotides are represented
        for nucleotide in 'ATCGN':
            if nucleotide not in composition:
                composition[nucleotide] = 0
        
        return dict(composition)
    
    def calculate_sequence_statistics(self, sequences: List) -> Dict[str, Any]:
        """
        Calculate comprehensive statistics for a list of sequences
        
        Args:
            sequences (List): List of Bio.SeqRecord objects
            
        Returns:
            Dict[str, Any]: Dictionary with sequence statistics
        """
        if not sequences:
            return {}
        
        lengths = [len(seq.seq) for seq in sequences]
        gc_contents = [self.calculate_gc_content(str(seq.seq)) for seq in sequences]
        
        stats = {
            'sequence_count': len(sequences),
            'total_length': sum(lengths),
            'average_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'length_std': np.std(lengths),
            'average_gc_content': np.mean(gc_contents),
            'median_gc_content': np.median(gc_contents),
            'min_gc_content': min(gc_contents),
            'max_gc_content': max(gc_contents),
            'gc_std': np.std(gc_contents)
        }
        
        return stats
    
    def calculate_quality_statistics(self, sequences: List) -> Optional[Dict[str, float]]:
        """
        Calculate quality score statistics for FASTQ sequences
        
        Args:
            sequences (List): List of Bio.SeqRecord objects with quality scores
            
        Returns:
            Optional[Dict[str, float]]: Quality statistics or None if no quality data
        """
        try:
            quality_scores = []
            for seq in sequences:
                if hasattr(seq, 'letter_annotations') and 'phred_quality' in seq.letter_annotations:
                    quality_scores.extend(seq.letter_annotations['phred_quality'])
            
            if not quality_scores:
                return None
            
            return {
                'average_quality': np.mean(quality_scores),
                'median_quality': np.median(quality_scores),
                'min_quality': min(quality_scores),
                'max_quality': max(quality_scores),
                'quality_std': np.std(quality_scores),
                'q20_bases': sum(1 for q in quality_scores if q >= 20),
                'q30_bases': sum(1 for q in quality_scores if q >= 30),
                'total_bases': len(quality_scores)
            }
        except Exception:
            return None
    
    def get_quality_data(self, sequences: List) -> Optional[List[List[float]]]:
        """
        Extract quality data for visualization
        
        Args:
            sequences (List): List of Bio.SeqRecord objects with quality scores
            
        Returns:
            Optional[List[List[float]]]: Quality data per sequence or None
        """
        try:
            quality_data = []
            for seq in sequences:
                if hasattr(seq, 'letter_annotations') and 'phred_quality' in seq.letter_annotations:
                    quality_data.append(seq.letter_annotations['phred_quality'])
            
            return quality_data if quality_data else None
        except Exception:
            return None
    
    def search_pattern(self, sequences: List, pattern: str) -> List[Dict[str, Any]]:
        """
        Search for a specific pattern in sequences
        
        Args:
            sequences (List): List of Bio.SeqRecord objects
            pattern (str): Pattern to search for
            
        Returns:
            List[Dict[str, Any]]: List of matches with details
        """
        if not pattern:
            return []
        
        pattern = pattern.upper()
        matches = []
        
        for i, seq in enumerate(sequences):
            sequence_str = str(seq.seq).upper()
            
            # Find all occurrences
            for match in re.finditer(pattern, sequence_str):
                matches.append({
                    'sequence_index': i,
                    'sequence_id': seq.id if hasattr(seq, 'id') else f"Sequence_{i}",
                    'start_position': match.start() + 1,  # 1-based indexing
                    'end_position': match.end(),
                    'pattern': pattern,
                    'context': sequence_str[max(0, match.start()-10):match.end()+10]
                })
        
        return matches
    
    def calculate_dinucleotide_frequency(self, sequence: str) -> Dict[str, int]:
        """
        Calculate dinucleotide frequency
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            Dict[str, int]: Dinucleotide frequencies
        """
        sequence = sequence.upper()
        dinucleotides = {}
        
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            if all(n in self.valid_nucleotides for n in dinuc):
                dinucleotides[dinuc] = dinucleotides.get(dinuc, 0) + 1
        
        return dinucleotides
    
    def find_orfs(self, sequence: str, min_length: int = 100) -> List[Dict[str, Any]]:
        """
        Find Open Reading Frames (ORFs) in a sequence
        
        Args:
            sequence (str): DNA sequence
            min_length (int): Minimum ORF length in nucleotides
            
        Returns:
            List[Dict[str, Any]]: List of ORFs with details
        """
        sequence = sequence.upper()
        start_codons = ['ATG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        orfs = []
        
        # Check all three reading frames
        for frame in range(3):
            for i in range(frame, len(sequence) - 2, 3):
                codon = sequence[i:i+3]
                
                if codon in start_codons:
                    # Look for stop codon
                    for j in range(i + 3, len(sequence) - 2, 3):
                        stop_codon = sequence[j:j+3]
                        if stop_codon in stop_codons:
                            orf_length = j + 3 - i
                            if orf_length >= min_length:
                                orfs.append({
                                    'start': i + 1,  # 1-based indexing
                                    'end': j + 3,
                                    'length': orf_length,
                                    'frame': frame + 1,
                                    'sequence': sequence[i:j+3]
                                })
                            break
        
        return orfs
    
    def reverse_complement(self, sequence: str) -> str:
        """
        Calculate reverse complement of a DNA sequence
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            str: Reverse complement sequence
        """
        sequence = sequence.upper()
        complement = ''.join(self.complement_map.get(base, base) for base in sequence)
        return complement[::-1]
    
    def calculate_complexity(self, sequence: str, window_size: int = 10) -> float:
        """
        Calculate sequence complexity using Shannon entropy
        
        Args:
            sequence (str): DNA sequence
            window_size (int): Window size for complexity calculation
            
        Returns:
            float: Average complexity score
        """
        if len(sequence) < window_size:
            return 0.0
        
        complexities = []
        
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]
            composition = Counter(window.upper())
            total = len(window)
            
            # Calculate Shannon entropy
            entropy = 0
            for count in composition.values():
                if count > 0:
                    p = count / total
                    entropy -= p * np.log2(p)
            
            complexities.append(entropy)
        
        return np.mean(complexities) if complexities else 0.0
