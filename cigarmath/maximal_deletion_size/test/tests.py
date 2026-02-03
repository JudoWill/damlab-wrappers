import os
import pytest # type: ignore
import yaml
import csv

def test_output_exists():
    """Test that output files were created"""
    assert os.path.exists('test_output/test_output.yaml')
    assert os.path.exists('test_output/test_read_level.csv')

def test_yaml_output_format():
    """Test that output YAML has correct structure and values"""
    with open('test_output/test_output.yaml') as f:
        metrics = yaml.safe_load(f)
    
    # Check required keys
    required_keys = {
        'sample_name',
        'mapped_reads',
        'average_ref_block_size',
        'average_largest_deletion_size',
        'total_reads'
    }
    assert set(metrics.keys()) == required_keys
    
    # Check sample name
    assert metrics['sample_name'] == 'test'
    
    # Check that metrics are non-negative
    assert metrics['mapped_reads'] >= 0
    assert metrics['average_ref_block_size'] >= 0
    assert metrics['average_largest_deletion_size'] >= 0
    assert metrics['total_reads'] >= 0

def test_csv_output_format():
    """Test that CSV output has correct structure and values"""
    with open('test_output/test_read_level.csv', 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    # Check CSV structure
    required_columns = {
        'read_name',
        'reference_block_start',
        'reference_block_stop',
        'reference_block_size',
        'largest_deletion_start',
        'largest_deletion_stop',
        'largest_deletion_size'
    }
    assert set(rows[0].keys()) == required_columns
    
    # Check data types and values
    for row in rows:
        # Check numeric fields
        assert int(row['reference_block_start']) >= 0
        assert int(row['reference_block_stop']) >= int(row['reference_block_start'])
        assert int(row['reference_block_size']) == int(row['reference_block_stop']) - int(row['reference_block_start'])
        
        # Check deletion fields
        if row['largest_deletion_start'] and row['largest_deletion_stop']:
            assert int(row['largest_deletion_start']) >= 0
            assert int(row['largest_deletion_stop']) >= int(row['largest_deletion_start'])
            assert int(row['largest_deletion_size']) == int(row['largest_deletion_stop']) - int(row['largest_deletion_start'])
        else:
            assert row['largest_deletion_size'] == '0'

def test_metrics_consistency():
    """Test that YAML metrics are consistent with CSV data"""
    with open('test_output/test_output.yaml') as f:
        metrics = yaml.safe_load(f)
    
    with open('test_output/test_read_level.csv', 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    # Count mapped reads
    mapped_reads = sum(1 for row in rows if int(row['reference_block_size']) > 0)
    assert mapped_reads == metrics['mapped_reads']
    
    # Calculate average reference block size
    if mapped_reads > 0:
        total_ref_size = sum(int(row['reference_block_size']) for row in rows)
        avg_ref_size = total_ref_size / mapped_reads
        assert abs(avg_ref_size - metrics['average_ref_block_size']) < 0.001
        
        # Calculate average largest deletion size
        total_del_size = sum(int(row['largest_deletion_size']) for row in rows)
        avg_del_size = total_del_size / mapped_reads
        assert abs(avg_del_size - metrics['average_largest_deletion_size']) < 0.001 