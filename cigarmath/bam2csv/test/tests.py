"""Unit tests for bam2csv functions"""

import csv

def test_output_format():
    """Test that output CSV has correct format and headers"""
    with open('test_output/test_output.csv', newline='') as f:
        reader = csv.reader(f)
        
        # Check header
        header = next(reader)
        expected_header = ["query_name", "reference_name", "reference_start", "query_sequence", "NM"]
        assert header == expected_header
        
        # Check data rows
        for row in reader:
            # Check number of fields matches header
            assert len(row) == len(expected_header)
            
            # Check types
            assert isinstance(row[0], str)  # query_name
            assert isinstance(row[1], str)  # reference_name
            assert row[2].isdigit()  # reference_start
            assert isinstance(row[3], str)  # query_sequence
            assert row[4].isdigit() or row[4] == ''  # NM tag

def test_row_count():
    """Test that the number of rows in CSV matches the number of records in BAM"""

    bam_count = 248
    
    # Count rows in CSV (excluding header)
    with open('test_output/test_output.csv', newline='') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        csv_count = sum(1 for _ in reader)
    
    # Check counts match
    assert bam_count == csv_count, f"BAM has {bam_count} records but CSV has {csv_count} rows" 