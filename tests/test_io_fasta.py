"""
Unit tests for the FASTA/FASTQ I/O module.

Run with: pytest tests/test_io_fasta.py -v
"""

import pytest
from app.core.io_fasta import (
    detect_format,
    parse_fasta,
    parse_fastq,
    parse_sequences,
    export_to_fasta,
    clean_sequence
)


class TestFormatDetection:
    """Test suite for file format detection."""
    
    def test_detect_fasta_by_content(self):
        """Test FASTA detection from content."""
        content = ">seq1\nATGC"
        assert detect_format(content) == 'fasta'
    
    def test_detect_fastq_by_content(self):
        """Test FASTQ detection from content."""
        content = "@seq1\nATGC\n+\n!!!!"
        assert detect_format(content) == 'fastq'
    
    def test_detect_fasta_by_extension(self):
        """Test FASTA detection from filename."""
        assert detect_format("content", "test.fasta") == 'fasta'
        assert detect_format("content", "test.fa") == 'fasta'
        assert detect_format("content", "test.fna") == 'fasta'
    
    def test_detect_fastq_by_extension(self):
        """Test FASTQ detection from filename."""
        assert detect_format("content", "test.fastq") == 'fastq'
        assert detect_format("content", "test.fq") == 'fastq'
    
    def test_detect_unknown_format(self):
        """Test handling of unknown format."""
        assert detect_format("random content") is None


class TestFastaParsing:
    """Test suite for FASTA parsing."""
    
    def test_parse_simple_fasta(self):
        """Test parsing simple FASTA."""
        content = ">seq1\nATGC\n>seq2\nGCTA"
        records = parse_fasta(content)
        
        assert len(records) == 2
        assert str(records[0].seq) == "ATGC"
        assert str(records[1].seq) == "GCTA"
    
    def test_parse_multiline_fasta(self):
        """Test parsing multi-line FASTA sequence."""
        content = ">seq1 description\nATGC\nGCTA\nAAAA"
        records = parse_fasta(content)
        
        assert len(records) == 1
        assert str(records[0].seq) == "ATGCGCTAAAA"
    
    def test_parse_fasta_with_description(self):
        """Test parsing FASTA with descriptions."""
        content = ">seq1 This is a description\nATGC"
        records = parse_fasta(content)
        
        assert records[0].id == "seq1"
        assert "description" in records[0].description
    
    def test_parse_invalid_fasta(self):
        """Test that invalid FASTA raises error."""
        content = "Not a FASTA file"
        
        with pytest.raises(ValueError):
            parse_fasta(content)


class TestFastqParsing:
    """Test suite for FASTQ parsing."""
    
    def test_parse_simple_fastq(self):
        """Test parsing simple FASTQ."""
        content = "@seq1\nATGC\n+\n!!!!"
        records = parse_fastq(content)
        
        assert len(records) == 1
        assert str(records[0].seq) == "ATGC"
    
    def test_parse_fastq_quality_scores(self):
        """Test that quality scores are parsed."""
        content = "@seq1\nATGC\n+\nIIII"
        records = parse_fastq(content)
        
        assert 'phred_quality' in records[0].letter_annotations
    
    def test_parse_invalid_fastq(self):
        """Test that invalid FASTQ raises error."""
        content = "Not a FASTQ file"
        
        with pytest.raises(ValueError):
            parse_fastq(content)


class TestExport:
    """Test suite for FASTA export."""
    
    def test_export_simple(self):
        """Test exporting sequences to FASTA."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        
        records = [
            SeqRecord(Seq("ATGC"), id="seq1"),
            SeqRecord(Seq("GCTA"), id="seq2")
        ]
        
        output = export_to_fasta(records)
        
        assert ">seq1" in output
        assert "ATGC" in output
        assert ">seq2" in output


class TestCleanSequence:
    """Test suite for sequence cleaning."""
    
    def test_clean_whitespace(self):
        """Test removal of whitespace."""
        sequence = "AT GC\nTA\tGC"
        assert clean_sequence(sequence) == "ATGCTAGC"
    
    def test_clean_uppercase(self):
        """Test conversion to uppercase."""
        sequence = "atgc"
        assert clean_sequence(sequence) == "ATGC"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
