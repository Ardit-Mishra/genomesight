from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
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
        Calculate comprehensive statistics for a list of sequences with advanced metrics
        
        Args:
            sequences (List): List of Bio.SeqRecord objects
            
        Returns:
            Dict[str, Any]: Dictionary with sequence statistics
        """
        if not sequences:
            return {}
        
        lengths = [len(seq.seq) for seq in sequences]
        gc_contents = [self.calculate_gc_content(str(seq.seq)) for seq in sequences]
        
        # Calculate additional advanced metrics
        complexities = [self.calculate_complexity(str(seq.seq)) for seq in sequences]
        
        # Nucleotide composition analysis
        total_composition = Counter()
        for seq in sequences:
            seq_comp = self.calculate_nucleotide_composition(str(seq.seq))
            for nucleotide, count in seq_comp.items():
                total_composition[nucleotide] += count
        
        total_bases = sum(total_composition.values())
        composition_percentages = {nuc: (count/total_bases)*100 for nuc, count in total_composition.items()}
        
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
            'gc_std': np.std(gc_contents),
            'average_complexity': np.mean(complexities),
            'total_nucleotide_composition': total_composition,
            'composition_percentages': composition_percentages,
            'at_content': composition_percentages.get('A', 0) + composition_percentages.get('T', 0),
            'purine_content': composition_percentages.get('A', 0) + composition_percentages.get('G', 0),
            'pyrimidine_content': composition_percentages.get('C', 0) + composition_percentages.get('T', 0)
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
                                    'sequence': sequence[i:j+3],
                                    'gc_content': self.calculate_gc_content(sequence[i:j+3])
                                })
                            break
        
        return orfs
    
    def analyze_all_advanced(self, sequences: List) -> Dict[str, Any]:
        """
        Run comprehensive advanced analysis on all sequences for research
        
        Args:
            sequences: List of Bio.SeqRecord objects
            
        Returns:
            Dict: Complete analysis results with research features
        """
        if not sequences:
            return {}
        
        results = {
            'basic_stats': self.calculate_sequence_statistics(sequences),
            'quality_stats': self.calculate_quality_statistics(sequences),
            'quality_data': self.get_quality_data(sequences)
        }
        
        # Add advanced research features
        try:
            results['motif_analysis'] = self.find_motifs(sequences)
        except Exception as e:
            results['motif_analysis'] = {'error': f'Motif analysis failed: {str(e)}'}
        
        try:
            results['codon_usage'] = self.analyze_codon_usage(sequences)
        except Exception as e:
            results['codon_usage'] = {'error': f'Codon analysis failed: {str(e)}'}
        
        try:
            results['orf_prediction'] = self.predict_orfs_advanced(sequences)
        except Exception as e:
            results['orf_prediction'] = {'error': f'ORF prediction failed: {str(e)}'}
        
        try:
            results['evolutionary_distances'] = self.calculate_evolutionary_distances(sequences)
        except Exception as e:
            results['evolutionary_distances'] = {'error': f'Distance calculation failed: {str(e)}'}
        
        try:
            results['consensus_analysis'] = self.perform_consensus_analysis(sequences)
        except Exception as e:
            results['consensus_analysis'] = {'error': f'Consensus analysis failed: {str(e)}'}
        
        return results
    
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
    
    def analyze_kmers(self, sequences: List, k: int = 3) -> Dict[str, Any]:
        """
        Analyze k-mer frequencies in sequences
        
        Args:
            sequences: List of Bio.SeqRecord objects
            k: K-mer length (default 3 for codons)
            
        Returns:
            Dict: K-mer analysis results
        """
        from collections import Counter
        import numpy as np
        
        all_kmers = Counter()
        kmer_stats = []
        
        for seq_record in sequences:
            sequence = str(seq_record.seq).upper()
            sequence_kmers = Counter()
            
            # Extract k-mers
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                if all(base in 'ATCG' for base in kmer):  # Only valid nucleotides
                    sequence_kmers[kmer] += 1
                    all_kmers[kmer] += 1
            
            # Calculate statistics for this sequence
            if sequence_kmers:
                total_kmers = sum(sequence_kmers.values())
                unique_kmers = len(sequence_kmers)
                most_common = sequence_kmers.most_common(10)
                
                # Calculate diversity metrics
                frequencies = np.array(list(sequence_kmers.values()))
                shannon_entropy = -np.sum((frequencies / total_kmers) * np.log2(frequencies / total_kmers))
                
                kmer_stats.append({
                    'sequence_id': seq_record.id,
                    'total_kmers': total_kmers,
                    'unique_kmers': unique_kmers,
                    'diversity': shannon_entropy,
                    'most_common': most_common,
                    'kmer_counts': sequence_kmers
                })
        
        # Global statistics
        total_global_kmers = sum(all_kmers.values())
        unique_global_kmers = len(all_kmers)
        global_most_common = all_kmers.most_common(20)
        
        # Calculate expected vs observed for GC bias
        gc_bias_analysis = self._analyze_gc_bias_in_kmers(all_kmers, k)
        
        return {
            'k': k,
            'global_stats': {
                'total_kmers': total_global_kmers,
                'unique_kmers': unique_global_kmers,
                'most_common': global_most_common,
                'all_kmers': all_kmers,
                'gc_bias': gc_bias_analysis
            },
            'per_sequence': kmer_stats
        }
    
    def _analyze_gc_bias_in_kmers(self, kmer_counts, k):
        """Analyze GC bias in k-mer distribution"""
        gc_rich_kmers = 0
        at_rich_kmers = 0
        balanced_kmers = 0
        
        for kmer, count in kmer_counts.items():
            gc_count = kmer.count('G') + kmer.count('C')
            gc_content = gc_count / len(kmer)
            
            if gc_content > 0.6:
                gc_rich_kmers += count
            elif gc_content < 0.4:
                at_rich_kmers += count
            else:
                balanced_kmers += count
        
        total = gc_rich_kmers + at_rich_kmers + balanced_kmers
        
        return {
            'gc_rich_percent': (gc_rich_kmers / total * 100) if total > 0 else 0,
            'at_rich_percent': (at_rich_kmers / total * 100) if total > 0 else 0,
            'balanced_percent': (balanced_kmers / total * 100) if total > 0 else 0
        }
    
    def analyze_reading_frames(self, sequences: List) -> Dict[str, Any]:
        """
        Comprehensive reading frame analysis
        
        Args:
            sequences: List of Bio.SeqRecord objects
            
        Returns:
            Dict: Reading frame analysis results
        """
        reading_frame_stats = []
        
        for seq_record in sequences:
            sequence = str(seq_record.seq).upper()
            seq_length = len(sequence)
            
            frame_analysis = {
                'sequence_id': seq_record.id,
                'sequence_length': seq_length,
                'frames': {}
            }
            
            # Analyze all 6 reading frames (3 forward, 3 reverse)
            for strand in ['+', '-']:
                if strand == '+':
                    seq_to_analyze = sequence
                else:
                    try:
                        seq_to_analyze = str(seq_record.seq.reverse_complement()).upper()
                    except Exception:
                        continue  # Skip if reverse complement not supported
                
                for frame in range(3):
                    frame_key = f"{strand}{frame + 1}"
                    
                    # Extract codons for this frame
                    codons = []
                    amino_acids = []
                    
                    for i in range(frame, len(seq_to_analyze) - 2, 3):
                        codon = seq_to_analyze[i:i+3]
                        if len(codon) == 3 and all(base in 'ATCG' for base in codon):
                            codons.append(codon)
                            # Simple codon to amino acid translation
                            aa = self._translate_codon(codon)
                            amino_acids.append(aa)
                    
                    # Analyze this reading frame
                    start_codons = sum(1 for codon in codons if codon == 'ATG')
                    stop_codons = sum(1 for codon in codons if codon in ['TAA', 'TAG', 'TGA'])
                    
                    # Calculate longest ORF in this frame
                    longest_orf = self._find_longest_orf_in_frame(codons)
                    
                    frame_analysis['frames'][frame_key] = {
                        'total_codons': len(codons),
                        'start_codons': start_codons,
                        'stop_codons': stop_codons,
                        'longest_orf': longest_orf,
                        'amino_acid_count': len([aa for aa in amino_acids if aa != '*']),
                        'stop_frequency': stop_codons / len(codons) if codons else 0
                    }
            
            reading_frame_stats.append(frame_analysis)
        
        return {
            'analysis_type': 'reading_frames',
            'total_sequences': len(sequences),
            'frame_stats': reading_frame_stats
        }
    
    def _translate_codon(self, codon):
        """Simple codon to amino acid translation"""
        codon_table = {
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
        return codon_table.get(codon, 'X')  # X for unknown
    
    def _find_longest_orf_in_frame(self, codons):
        """Find the longest ORF in a list of codons"""
        longest = 0
        current_length = 0
        in_orf = False
        
        for codon in codons:
            if codon == 'ATG':  # Start codon
                if not in_orf:
                    in_orf = True
                    current_length = 3  # Length of start codon
            elif codon in ['TAA', 'TAG', 'TGA']:  # Stop codons
                if in_orf:
                    current_length += 3  # Add stop codon length
                    longest = max(longest, current_length)
                    in_orf = False
                    current_length = 0
            elif in_orf:
                current_length += 3
        
        return longest

    def find_motifs(self, sequences: List, min_length: int = 6, max_length: int = 12) -> Dict[str, Any]:
        """
        Discover conserved motifs in sequences using k-mer analysis
        
        Args:
            sequences: List of Bio.SeqRecord objects
            min_length: Minimum motif length
            max_length: Maximum motif length
            
        Returns:
            Dict: Motif analysis results
        """
        motifs = {}
        
        for k in range(min_length, max_length + 1):
            kmer_counts = Counter()
            total_kmers = 0
            
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                for i in range(len(seq_str) - k + 1):
                    kmer = seq_str[i:i+k]
                    if all(base in 'ATCG' for base in kmer):
                        kmer_counts[kmer] += 1
                        total_kmers += 1
            
            # Find enriched motifs (appearing in multiple sequences)
            enriched_motifs = []
            for kmer, count in kmer_counts.most_common(20):
                # Count how many sequences contain this motif
                seq_count = sum(1 for seq in sequences if kmer in str(seq.seq).upper())
                frequency = count / total_kmers if total_kmers > 0 else 0
                conservation = seq_count / len(sequences) if sequences else 0
                
                if conservation >= 0.3:  # Present in at least 30% of sequences
                    enriched_motifs.append({
                        'motif': kmer,
                        'count': count,
                        'frequency': frequency,
                        'conservation': conservation,
                        'sequences_with_motif': seq_count
                    })
            
            if enriched_motifs:
                motifs[f'{k}-mer'] = enriched_motifs[:10]  # Top 10 motifs
        
        return motifs
    
    def analyze_codon_usage(self, sequences: List) -> Dict[str, Any]:
        """
        Analyze codon usage bias in coding sequences
        
        Args:
            sequences: List of Bio.SeqRecord objects
            
        Returns:
            Dict: Codon usage analysis results
        """
        codon_usage = Counter()
        amino_acid_usage = Counter()
        total_codons = 0
        
        # Standard genetic code
        genetic_code = {
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
        
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            # Analyze in reading frames
            for frame in range(3):
                for i in range(frame, len(seq_str) - 2, 3):
                    codon = seq_str[i:i+3]
                    if len(codon) == 3 and all(base in 'ATCG' for base in codon):
                        codon_usage[codon] += 1
                        total_codons += 1
                        if codon in genetic_code:
                            amino_acid = genetic_code[codon]
                            amino_acid_usage[amino_acid] += 1
        
        # Calculate relative synonymous codon usage (RSCU)
        rscu_values = {}
        for aa in set(genetic_code.values()):
            if aa != '*':  # Skip stop codons
                aa_codons = [codon for codon, amino in genetic_code.items() if amino == aa]
                aa_total = sum(codon_usage[codon] for codon in aa_codons)
                if aa_total > 0:
                    for codon in aa_codons:
                        expected = aa_total / len(aa_codons)
                        rscu = (codon_usage[codon] / expected) if expected > 0 else 0
                        rscu_values[codon] = rscu
        
        return {
            'total_codons': total_codons,
            'codon_frequencies': dict(codon_usage.most_common()),
            'amino_acid_frequencies': dict(amino_acid_usage.most_common()),
            'rscu_values': rscu_values,
            'most_biased_codons': sorted(rscu_values.items(), key=lambda x: abs(x[1] - 1), reverse=True)[:10]
        }
    
    def predict_orfs_advanced(self, sequences: List, min_length: int = 300) -> Dict[str, Any]:
        """
        Advanced ORF prediction with protein analysis
        
        Args:
            sequences: List of Bio.SeqRecord objects
            min_length: Minimum ORF length in base pairs
            
        Returns:
            Dict: ORF prediction results with protein analysis
        """
        orfs_found = []
        start_codons = ['ATG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        for seq_idx, seq in enumerate(sequences):
            seq_str = str(seq.seq).upper()
            
            # Check all 6 reading frames (3 forward, 3 reverse)
            for strand in [1, -1]:
                current_seq = seq_str
                if strand == -1:
                    current_seq = str(seq.seq.reverse_complement()).upper()
                
                for frame in range(3):
                    i = frame
                    while i < len(current_seq) - 2:
                        codon = current_seq[i:i+3]
                        
                        # Look for start codon
                        if codon in start_codons:
                            start_pos = i
                            j = i + 3
                            
                            # Look for stop codon
                            while j < len(current_seq) - 2:
                                stop_codon = current_seq[j:j+3]
                                if stop_codon in stop_codons:
                                    orf_length = j - start_pos + 3
                                    if orf_length >= min_length:
                                        orf_seq = current_seq[start_pos:j+3]
                                        
                                        # Translate to protein
                                        protein_seq = str(Seq(orf_seq).translate())
                                        
                                        # Analyze protein properties
                                        try:
                                            protein_analysis = ProteinAnalysis(protein_seq[:-1])  # Remove stop codon
                                            molecular_weight = protein_analysis.molecular_weight()
                                            isoelectric_point = protein_analysis.isoelectric_point()
                                        except:
                                            molecular_weight = 0
                                            isoelectric_point = 0
                                        
                                        orfs_found.append({
                                            'sequence_id': seq_idx,
                                            'sequence_name': seq.id if hasattr(seq, 'id') else f'seq_{seq_idx}',
                                            'strand': '+' if strand == 1 else '-',
                                            'frame': frame + 1,
                                            'start': start_pos + 1,
                                            'end': j + 3,
                                            'length': orf_length,
                                            'nucleotide_sequence': orf_seq,
                                            'protein_sequence': protein_seq,
                                            'protein_length': len(protein_seq) - 1,
                                            'molecular_weight': molecular_weight,
                                            'isoelectric_point': isoelectric_point
                                        })
                                    
                                    i = j + 3
                                    break
                                j += 3
                            else:
                                i += 3
                        else:
                            i += 3
        
        return {
            'total_orfs': len(orfs_found),
            'orfs': orfs_found,
            'avg_orf_length': np.mean([orf['length'] for orf in orfs_found]) if orfs_found else 0,
            'longest_orf': max(orfs_found, key=lambda x: x['length']) if orfs_found else None,
            'avg_protein_length': np.mean([orf['protein_length'] for orf in orfs_found]) if orfs_found else 0
        }
    
    def calculate_evolutionary_distances(self, sequences: List) -> Dict[str, Any]:
        """
        Calculate pairwise evolutionary distances between sequences
        
        Args:
            sequences: List of Bio.SeqRecord objects
            
        Returns:
            Dict: Distance matrix and phylogenetic analysis
        """
        if len(sequences) < 2:
            return {'error': 'Need at least 2 sequences for distance calculation'}
        
        # Limit to reasonable number for performance
        sample_sequences = sequences[:50] if len(sequences) > 50 else sequences
        n_seqs = len(sample_sequences)
        
        distance_matrix = np.zeros((n_seqs, n_seqs))
        sequence_names = []
        
        for i, seq in enumerate(sample_sequences):
            sequence_names.append(seq.id if hasattr(seq, 'id') else f'seq_{i}')
        
        # Calculate pairwise distances using p-distance
        for i in range(n_seqs):
            for j in range(i + 1, n_seqs):
                seq1 = str(sample_sequences[i].seq).upper()
                seq2 = str(sample_sequences[j].seq).upper()
                
                # Align sequences to same length for comparison
                min_len = min(len(seq1), len(seq2))
                if min_len > 0:
                    seq1_trimmed = seq1[:min_len]
                    seq2_trimmed = seq2[:min_len]
                    
                    # Calculate p-distance (proportion of differing sites)
                    differences = sum(1 for a, b in zip(seq1_trimmed, seq2_trimmed) if a != b)
                    distance = differences / min_len
                    
                    distance_matrix[i][j] = distance
                    distance_matrix[j][i] = distance
        
        # Find most similar and divergent pairs
        non_zero_distances = distance_matrix[distance_matrix > 0]
        most_similar_pairs = []
        most_divergent_pairs = []
        
        if len(non_zero_distances) > 0:
            min_dist = np.min(non_zero_distances)
            max_dist = np.max(distance_matrix)
            
            # Find pairs with minimum and maximum distances
            for i in range(n_seqs):
                for j in range(i + 1, n_seqs):
                    if distance_matrix[i][j] == min_dist:
                        most_similar_pairs.append((sequence_names[i], sequence_names[j], min_dist))
                    if distance_matrix[i][j] == max_dist:
                        most_divergent_pairs.append((sequence_names[i], sequence_names[j], max_dist))
        
        return {
            'distance_matrix': distance_matrix.tolist(),
            'sequence_names': sequence_names,
            'avg_distance': np.mean(non_zero_distances) if len(non_zero_distances) > 0 else 0,
            'max_distance': np.max(distance_matrix),
            'min_distance': np.min(non_zero_distances) if len(non_zero_distances) > 0 else 0,
            'most_similar_pairs': most_similar_pairs[:5],
            'most_divergent_pairs': most_divergent_pairs[:5]
        }
    
    def perform_consensus_analysis(self, sequences: List) -> Dict[str, Any]:
        """
        Generate consensus sequence and variability analysis
        
        Args:
            sequences: List of Bio.SeqRecord objects
            
        Returns:
            Dict: Consensus sequence and variability data
        """
        if not sequences:
            return {'error': 'No sequences provided'}
        
        # Get sequences as strings
        seq_strings = [str(seq.seq).upper() for seq in sequences]
        max_length = max(len(s) for s in seq_strings)
        
        # Pad sequences to same length
        padded_sequences = [s.ljust(max_length, 'N') for s in seq_strings]
        
        consensus = []
        variability = []
        
        for pos in range(max_length):
            # Count nucleotides at this position
            nucleotides = [seq[pos] for seq in padded_sequences if pos < len(seq) and seq[pos] != 'N']
            
            if nucleotides:
                counts = Counter(nucleotides)
                most_common = counts.most_common(1)[0]
                consensus_base = most_common[0]
                consensus_freq = most_common[1] / len(nucleotides)
                
                consensus.append(consensus_base)
                variability.append(1 - consensus_freq)  # Variability = 1 - conservation
            else:
                consensus.append('N')
                variability.append(0)
        
        consensus_sequence = ''.join(consensus)
        avg_variability = np.mean(variability)
        
        # Find most variable positions
        variable_positions = [(i+1, variability[i]) for i in range(len(variability))]
        variable_positions.sort(key=lambda x: x[1], reverse=True)
        
        return {
            'consensus_sequence': consensus_sequence,
            'consensus_length': len(consensus_sequence),
            'average_variability': avg_variability,
            'most_variable_positions': variable_positions[:20],
            'conservation_profile': variability
        }
