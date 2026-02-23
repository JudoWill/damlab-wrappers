import os
import yaml
import csv


def test_output_files_exist():
    """Test that all output files were created"""
    assert os.path.exists('test_output/reads.csv'), "reads.csv not found"
    assert os.path.exists('test_output/deletions.csv'), "deletions.csv not found"
    assert os.path.exists('test_output/summary.yaml'), "summary.yaml not found"


def test_reads_csv_structure():
    """Test that reads CSV has correct column structure"""
    with open('test_output/reads.csv', 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    expected_columns = {'read_name', 'reference_start', 'reference_end', 'deletions'}
    assert set(rows[0].keys()) == expected_columns, f"Expected columns {expected_columns}, got {set(rows[0].keys())}"
    
    assert len(rows) > 0, "reads.csv should have at least one row"


def test_deletions_csv_structure():
    """Test that deletions CSV has correct column structure"""
    with open('test_output/deletions.csv', 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    expected_columns = {'deletion_start', 'deletion_end', 'deletion_size', 'read_count', 'coverage_count'}
    
    if len(rows) > 0:
        assert set(rows[0].keys()) == expected_columns, f"Expected columns {expected_columns}, got {set(rows[0].keys())}"
        
        for row in rows:
            assert int(row['deletion_size']) == int(row['deletion_end']) - int(row['deletion_start']), \
                "deletion_size should equal deletion_end - deletion_start"
            assert int(row['read_count']) > 0, "read_count should be positive"
            assert int(row['coverage_count']) >= int(row['read_count']), \
                "coverage_count should be >= read_count"


def test_summary_yaml_structure():
    """Test that summary YAML has correct structure"""
    with open('test_output/summary.yaml') as f:
        content = f.read()
        assert content.startswith('# Cigarmath Deletion Block Detection'), \
            "YAML should start with MultiQC marker comment"
        
        summary = yaml.safe_load(content)
    
    required_keys = {
        'sample_name',
        'total_reads',
        'reads_with_deletions',
        'unique_deletion_count',
        'total_deletion_count',
        'deletion_frequency',
        'min_deletion_size',
        'input_bam_count',
        'allowedlist_used',
        'allowedlist_size'
    }
    
    assert set(summary.keys()) == required_keys, \
        f"Expected keys {required_keys}, got {set(summary.keys())}"


def test_summary_values():
    """Test that summary values are reasonable"""
    with open('test_output/summary.yaml') as f:
        summary = yaml.safe_load(f)
    
    assert summary['sample_name'] == 'test', "sample_name should be 'test'"
    assert summary['total_reads'] > 0, "total_reads should be positive"
    assert summary['reads_with_deletions'] >= 0, "reads_with_deletions should be non-negative"
    assert summary['reads_with_deletions'] <= summary['total_reads'], \
        "reads_with_deletions should not exceed total_reads"
    assert summary['unique_deletion_count'] >= 0, "unique_deletion_count should be non-negative"
    assert summary['total_deletion_count'] >= 0, "total_deletion_count should be non-negative"
    assert 0.0 <= summary['deletion_frequency'] <= 1.0, \
        "deletion_frequency should be between 0 and 1"
    assert summary['min_deletion_size'] == 50, "min_deletion_size should be 50"
    assert summary['input_bam_count'] == 1, "input_bam_count should be 1"
    assert summary['allowedlist_used'] == False, "allowedlist_used should be False"
    assert summary['allowedlist_size'] == 0, "allowedlist_size should be 0"


def test_cross_check_counts():
    """Cross-check counts between reads CSV and summary YAML"""
    with open('test_output/summary.yaml') as f:
        summary = yaml.safe_load(f)
    
    with open('test_output/reads.csv', 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    assert len(rows) == summary['total_reads'], \
        f"Number of rows in reads.csv ({len(rows)}) should match total_reads ({summary['total_reads']})"
    
    reads_with_dels = sum(1 for row in rows if row['deletions'])
    assert reads_with_dels == summary['reads_with_deletions'], \
        f"Reads with deletions ({reads_with_dels}) should match summary ({summary['reads_with_deletions']})"


def test_deletion_counts_match():
    """Verify deletion counts in deletions.csv sum to total_deletion_count in summary"""
    with open('test_output/summary.yaml') as f:
        summary = yaml.safe_load(f)
    
    with open('test_output/deletions.csv', 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    total_from_csv = sum(int(row['read_count']) for row in rows)
    assert total_from_csv == summary['total_deletion_count'], \
        f"Sum of read_count ({total_from_csv}) should match total_deletion_count ({summary['total_deletion_count']})"
    
    assert len(rows) == summary['unique_deletion_count'], \
        f"Number of deletion rows ({len(rows)}) should match unique_deletion_count ({summary['unique_deletion_count']})"
