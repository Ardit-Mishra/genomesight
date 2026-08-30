"""
Unit tests for the FASTA/FASTQ I/O module.

Run with: pytest tests/test_io_fasta.py -v
"""

import pytest
from app.core.io_fasta import (
    detect_format,
    parse_fasta,
    parse_fastq,
    parse_genbank,
    parse_sequences,
    export_to_fasta,
    extract_genbank_features,
    clean_sequence
)


GENBANK_SINGLE = """LOCUS       test_seq                  30 bp    DNA     linear   UNK 01-JAN-2024
DEFINITION  synthetic test record.
ACCESSION   test_seq
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="synthetic organism"
     CDS             1..30
                     /gene="testgene"
                     /product="test protein"
ORIGIN
        1 atgaaatagt gaatgcatgc atgcatgcat
//
"""

GENBANK_MULTI = """LOCUS       seq_one                   12 bp    DNA     linear   UNK 01-JAN-2024
DEFINITION  first record.
ACCESSION   seq_one
FEATURES             Location/Qualifiers
     source          1..12
                     /organism="synthetic organism"
ORIGIN
        1 atgaaatagt ga
//
LOCUS       seq_two                   12 bp    DNA     linear   UNK 01-JAN-2024
DEFINITION  second record, no features.
ACCESSION   seq_two
ORIGIN
        1 gcgcgcgcgc gc
//
"""


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
        # ATGC + GCTA + AAAA = 12 residues. The previous expectation dropped a
        # trailing A (11 residues), so this test failed against a correct parser.
        assert str(records[0].seq) == "ATGCGCTAAAAA"
        assert len(records[0].seq) == 12
    
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


class TestGenBankFormatDetection:
    """Test suite for GenBank format detection."""

    def test_detect_genbank_by_content(self):
        assert detect_format(GENBANK_SINGLE) == 'genbank'

    def test_detect_genbank_by_extension(self):
        assert detect_format("content", "test.gb") == 'genbank'
        assert detect_format("content", "test.gbk") == 'genbank'
        assert detect_format("content", "test.genbank") == 'genbank'

    def test_content_based_detection_beats_a_misleading_extension(self):
        # A LOCUS-first payload is GenBank even if someone names the file
        # wrong — content inspection runs before the extension fallback.
        assert detect_format(GENBANK_SINGLE, "misnamed.fasta") == 'genbank'


class TestGenBankParsing:
    """Test suite for GenBank parsing, including feature annotations."""

    def test_parse_single_record(self):
        records = parse_genbank(GENBANK_SINGLE)
        assert len(records) == 1
        assert str(records[0].seq) == "ATGAAATAGTGAATGCATGCATGCATGCAT"

    def test_parse_multi_record(self):
        records = parse_genbank(GENBANK_MULTI)
        assert len(records) == 2
        assert str(records[0].seq) == "ATGAAATAGTGA"
        assert str(records[1].seq) == "GCGCGCGCGCGC"

    def test_parse_invalid_genbank_raises(self):
        with pytest.raises(ValueError):
            parse_genbank("Not a GenBank file")

    def test_parse_sequences_routes_genbank_through_auto_detection(self):
        records, file_format = parse_sequences(GENBANK_SINGLE, "test.gb")
        assert file_format == 'genbank'
        assert len(records) == 1


class TestGenBankFeatureExtraction:
    """Feature annotations must survive parsing and flatten cleanly."""

    def test_features_are_extracted_with_qualifiers(self):
        records = parse_genbank(GENBANK_SINGLE)
        features = extract_genbank_features(records[0])

        assert len(features) == 2  # source + CDS
        cds = next(f for f in features if f['type'] == 'CDS')
        assert cds['start'] == 1
        assert cds['end'] == 30
        assert cds['strand'] == '+'
        assert cds['qualifiers']['gene'] == ['testgene']
        assert cds['qualifiers']['product'] == ['test protein']

    def test_record_with_no_features_returns_empty_list_not_an_error(self):
        records = parse_genbank(GENBANK_MULTI)
        second_record_features = extract_genbank_features(records[1])
        assert second_record_features == []

    def test_fasta_record_has_no_features(self):
        # FASTA records never carry Bio.SeqFeature annotations — confirm
        # this returns an empty list rather than raising.
        records = parse_fasta(">seq1\nATGC")
        assert extract_genbank_features(records[0]) == []


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
