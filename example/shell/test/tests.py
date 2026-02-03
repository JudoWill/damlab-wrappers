"""Unit tests for seqkit translate wrapper"""

import os
import pytest
from Bio import SeqIO

def test_output_exists():
    """Test that output files were created"""
    assert os.path.exists('test_output/basic.fasta'), "Basic output not found"
    assert os.path.exists('test_output/frame2.fasta'), "Frame 2 output not found"
    assert os.path.exists('test_output/table2.fasta'), "Table 2 output not found"

def test_basic_translation():
    """Test basic translation with default parameters"""
    records = list(SeqIO.parse('test_output/basic.fasta', 'fasta'))
    assert len(records) == 2, "Expected 2 sequences in output"
    
    # Check first sequence
    assert records[0].id == "seq1", "Expected sequence ID to be preserved"
    assert str(records[0].seq) == "MAMAK", "Incorrect translation"
    
    # Check second sequence
    assert records[1].id == "seq2", "Expected sequence ID to be preserved"
    assert str(records[1].seq) == "MAMAK", "Incorrect translation"

def test_frame1_translation():
    """Test translation with frame 1"""
    records = list(SeqIO.parse('test_output/frame2.fasta', 'fasta'))
    assert len(records) == 2, "Expected 2 sequences in output"
    
    # Check first sequence
    assert records[0].id == "seq1", "Expected sequence ID to be preserved"
    assert str(records[0].seq) == "WPWP", "Incorrect translation in frame 2"
    
    # Check second sequence
    assert records[1].id == "seq2", "Expected sequence ID to be preserved"
    assert str(records[1].seq) == "WPWPX", "Incorrect translation in frame 2"

def test_table2_translation():
    """Test translation with table 2 (vertebrate mitochondrial)"""
    records = list(SeqIO.parse('test_output/table2.fasta', 'fasta'))
    assert len(records) == 2, "Expected 2 sequences in output"
    
    # Check first sequence
    assert records[0].id == "seq1", "Expected sequence ID to be preserved"
    # In table 2, ATG->M, GCC->A, ATG->M, GCC->A, AAA->K
    assert str(records[0].seq) == "MAMAK", "Incorrect translation with table 2"
    
    # Check second sequence
    assert records[1].id == "seq2", "Expected sequence ID to be preserved"
    assert str(records[1].seq) == "MAMAK", "Incorrect translation with table 2" 