"""Unit tests for pileup functions"""

import os
import pytest
import csv
import math

def test_output_exists():
    """Test that output TSV file was created"""
    assert os.path.exists('test_output/test_output.tsv')

def test_output_format():
    """Test that output TSV has correct format"""
    with open('test_output/test_output.tsv', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        
        # Check header
        header = next(reader)
        expected_header = ['Position', 'Depth', 'Entropy', 'A', 'C', 'G', 'T', 'N', '-']
        assert header == expected_header
        
        # Check data rows
        for row in reader:
            # Check number of fields
            assert len(row) == 9
            
            # Check types
            assert row[0].isdigit()  # Position
            assert row[1].isdigit()  # Depth
            assert 0 <= float(row[2]) <= math.log2(6)  # Entropy (max for 6 possibilities)
            assert all(field.isdigit() for field in row[3:])  # Base counts

def test_depth_calculation():
    """Test that depths are calculated correctly"""
    with open('test_output/test_output.tsv', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # Skip header
        
        for row in reader:
            depth = int(row[1])
            # Sum all non-gap bases (A,C,G,T,N)
            calculated_depth = sum(int(count) for count in row[3:-1])
            assert depth == calculated_depth

def test_entropy_calculation():
    """Test that entropy values are reasonable"""
    with open('test_output/test_output.tsv', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # Skip header
        
        for row in reader:
            entropy = float(row[2])
            # Entropy should be between 0 and log2(6) ≈ 2.58
            assert 0 <= entropy <= math.log2(6)
            
            # If all counts are the same base, entropy should be 0
            counts = [int(x) for x in row[3:]]
            if sum(1 for x in counts if x > 0) == 1:
                assert entropy == 0 