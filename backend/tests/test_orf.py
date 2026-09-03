"""
Unit tests for the ORF detection module.

Run with: pytest tests/test_orf.py -v
"""

import pytest
from app.core.orf import (
    find_orfs,
    translate_sequence,
    reverse_complement,
    calculate_gc,
    get_orf_summary
)


class TestTranslation:
    """Test suite for sequence translation."""
    
    def test_translate_simple(self):
        """Test simple translation."""
        sequence = "ATGAAATAG"
        protein = translate_sequence(sequence)
        assert protein == "MK*"
    
    def test_translate_all_frames(self):
        """Test translation in different frames."""
        sequence = "AATGAAATAG"
        
        frame0 = translate_sequence(sequence, frame=0)
        frame1 = translate_sequence(sequence, frame=1)
        
        assert frame0 != frame1
    
    def test_translate_stop_codon(self):
        """Test that stop codons translate to *."""
        assert translate_sequence("TAA")[0] == "*"
        assert translate_sequence("TAG")[0] == "*"
        assert translate_sequence("TGA")[0] == "*"


class TestReverseComplement:
    """Test suite for reverse complement."""
    
    def test_reverse_complement_simple(self):
        """Test simple reverse complement."""
        assert reverse_complement("ATGC") == "GCAT"
    
    def test_reverse_complement_symmetry(self):
        """Test that double reverse complement returns original."""
        sequence = "ATGCATGC"
        assert reverse_complement(reverse_complement(sequence)) == sequence
    
    def test_reverse_complement_with_n(self):
        """Test reverse complement with N bases."""
        assert reverse_complement("ATNGC") == "GCNAT"


class TestGCCalculation:
    """Test suite for GC content calculation."""
    
    def test_gc_balanced(self):
        """Test balanced GC content."""
        assert calculate_gc("ATGC") == 50.0
    
    def test_gc_all_gc(self):
        """Test 100% GC content."""
        assert calculate_gc("GCGC") == 100.0
    
    def test_gc_no_gc(self):
        """Test 0% GC content."""
        assert calculate_gc("ATAT") == 0.0


class TestORFDetection:
    """Test suite for ORF finding."""
    
    def test_find_simple_orf(self):
        """Test finding a simple ORF."""
        sequence = "ATGAAAAAAAAATGA"
        orfs = find_orfs(sequence, min_length=9, include_reverse=False)
        
        assert len(orfs) == 1
        assert orfs[0].start_codon == "ATG"
        assert orfs[0].stop_codon == "TGA"
    
    def test_find_multiple_orfs(self):
        """Test finding multiple ORFs."""
        sequence = "ATGAAATAAATGGGGTGA"
        orfs = find_orfs(sequence, min_length=9, include_reverse=False)
        
        assert len(orfs) >= 1
    
    def test_find_no_orf(self):
        """Test when no ORF exists."""
        sequence = "AAAAAAAAAA"
        orfs = find_orfs(sequence, min_length=9)
        
        assert len(orfs) == 0
    
    def test_min_length_filter(self):
        """Test that minimum length filter works."""
        sequence = "ATGAAATAA"
        orfs_short = find_orfs(sequence, min_length=9, include_reverse=False)
        orfs_long = find_orfs(sequence, min_length=100, include_reverse=False)
        
        assert len(orfs_short) == 1
        assert len(orfs_long) == 0
    
    def test_reverse_complement_orfs(self):
        """Test finding ORFs on reverse strand."""
        sequence = "ATGAAATAA"
        orfs_forward = find_orfs(sequence, min_length=9, include_reverse=False)
        orfs_both = find_orfs(sequence, min_length=9, include_reverse=True)
        
        assert len(orfs_both) >= len(orfs_forward)
    
    def test_orf_attributes(self):
        """Test that ORF objects have expected attributes."""
        sequence = "ATGAAAAAAAAATGA"
        orfs = find_orfs(sequence, min_length=9, include_reverse=False)
        
        orf = orfs[0]
        assert hasattr(orf, 'start')
        assert hasattr(orf, 'end')
        assert hasattr(orf, 'length_nt')
        assert hasattr(orf, 'protein')
        assert hasattr(orf, 'gc_content')


class TestORFSummary:
    """Test suite for ORF summary statistics."""
    
    def test_summary_empty(self):
        """Test summary with no ORFs."""
        summary = get_orf_summary([])
        assert summary['total'] == 0
    
    def test_summary_with_orfs(self):
        """Test summary calculation."""
        sequence = "ATGAAAAAAAAATGA"
        orfs = find_orfs(sequence, min_length=9)
        summary = get_orf_summary(orfs)
        
        if orfs:
            assert summary['total'] == len(orfs)
            assert 'average_length' in summary
            assert 'max_length' in summary


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
