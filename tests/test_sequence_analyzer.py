"""
Unit tests for the SequenceAnalyzer module.

Run with: pytest tests/test_sequence_analyzer.py -v
"""

import pytest
from app.core.sequence_analyzer import SequenceAnalyzer


class TestSequenceAnalyzer:
    """Test suite for SequenceAnalyzer class."""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer instance for tests."""
        return SequenceAnalyzer()
    
    def test_gc_content_balanced(self, analyzer):
        """Test GC content with balanced sequence."""
        sequence = "ATGC"
        assert analyzer.calculate_gc_content(sequence) == 50.0
    
    def test_gc_content_all_gc(self, analyzer):
        """Test GC content with all GC sequence."""
        sequence = "GCGCGCGC"
        assert analyzer.calculate_gc_content(sequence) == 100.0
    
    def test_gc_content_no_gc(self, analyzer):
        """Test GC content with no GC bases."""
        sequence = "ATATAT"
        assert analyzer.calculate_gc_content(sequence) == 0.0
    
    def test_gc_content_empty(self, analyzer):
        """Test GC content with empty sequence."""
        assert analyzer.calculate_gc_content("") == 0.0
    
    def test_gc_content_case_insensitive(self, analyzer):
        """Test that GC calculation is case insensitive."""
        assert analyzer.calculate_gc_content("atgc") == 50.0
        assert analyzer.calculate_gc_content("ATGC") == 50.0
        assert analyzer.calculate_gc_content("AtGc") == 50.0
    
    def test_nucleotide_composition(self, analyzer):
        """Test nucleotide counting."""
        sequence = "AATTCCGG"
        composition = analyzer.calculate_nucleotide_composition(sequence)
        
        assert composition['A'] == 2
        assert composition['T'] == 2
        assert composition['C'] == 2
        assert composition['G'] == 2
    
    def test_nucleotide_composition_with_n(self, analyzer):
        """Test nucleotide counting with ambiguous bases."""
        sequence = "ATCGNN"
        composition = analyzer.calculate_nucleotide_composition(sequence)
        
        assert composition['N'] == 2
        assert composition['A'] == 1
    
    def test_reverse_complement(self, analyzer):
        """Test reverse complement generation."""
        sequence = "ATGC"
        assert analyzer.reverse_complement(sequence) == "GCAT"
    
    def test_reverse_complement_longer(self, analyzer):
        """Test reverse complement with longer sequence."""
        sequence = "AATTCCGG"
        assert analyzer.reverse_complement(sequence) == "CCGGAATT"
    
    def test_complexity_low(self, analyzer):
        """Test complexity for repetitive sequence."""
        sequence = "AAAAAAAAAA"
        complexity = analyzer.calculate_complexity(sequence, window_size=5)
        assert complexity < 1.0
    
    def test_complexity_high(self, analyzer):
        """Test complexity for diverse sequence."""
        sequence = "ATCGATCGATCGATCGATCG"
        complexity = analyzer.calculate_complexity(sequence, window_size=5)
        assert complexity > 1.0
    
    def test_translate_codon(self, analyzer):
        """Test single codon translation."""
        assert analyzer.translate_codon("ATG") == "M"
        assert analyzer.translate_codon("TAA") == "*"
        assert analyzer.translate_codon("TTT") == "F"
    
    def test_translate_codon_case_insensitive(self, analyzer):
        """Test that codon translation is case insensitive."""
        assert analyzer.translate_codon("atg") == "M"
        assert analyzer.translate_codon("Atg") == "M"


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    @pytest.fixture
    def analyzer(self):
        return SequenceAnalyzer()
    
    def test_composition_empty_sequence(self, analyzer):
        """Test composition with empty sequence."""
        composition = analyzer.calculate_nucleotide_composition("")
        assert composition['A'] == 0
        assert composition['T'] == 0
    
    def test_complexity_short_sequence(self, analyzer):
        """Test complexity with sequence shorter than window."""
        sequence = "ATG"
        complexity = analyzer.calculate_complexity(sequence, window_size=10)
        assert complexity == 0.0
    
    def test_unknown_codon(self, analyzer):
        """Test translation of unknown codon."""
        assert analyzer.translate_codon("XXX") == "X"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
