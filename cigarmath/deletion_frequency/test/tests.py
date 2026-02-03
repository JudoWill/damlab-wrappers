import os
import yaml
import csv

def test_output_exists():
    """Test that output YAML file was created"""
    assert os.path.exists('test_output/test_output.yaml')
    assert os.path.exists('test_output/test_read_level.csv')

def test_output_format():
    """Test that output YAML has correct structure"""
    with open('test_output/test_output.yaml') as f:
        stats = yaml.safe_load(f)
    
    required_keys = {
        'full_deletion_frequency',
        'partial_deletion_frequency',
        'reads_covered_full_deleted',
        'reads_covered_partial_deleted',
        'reads_covering_required',
        'region_name',
        'sample_name',
        'total_reads',
    }
    
    assert set(stats.keys()) == required_keys
    
    correct_values = {
        'full_deletion_frequency': 0.43333333333333335,
        'partial_deletion_frequency': 0.36666666666666664,
        'reads_covered_full_deleted': 13,
        'reads_covered_partial_deleted': 11,
        'reads_covering_required': 30,
        'region_name': 'region',
        'sample_name': 'test',
        'total_reads': 180
    }

    for key, value in correct_values.items():
        if isinstance(value, str):
            assert stats[key] == value, f'{key} is not equal to {value}'
        else:
            assert abs(stats[key] - value) < 0.001, f'{key} is not equal to {value}'


def test_csv_output_format():
    """Test that CSV output has correct structure and values"""
    with open('test_output/test_read_level.csv', 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    # Check CSV structure
    assert set(rows[0].keys()) == {'read_name','covers_required','full_deleted','partial_deleted'}
    
    # Count rows with covers_required=True and has_deletion=True
    covers_required_count = sum(1 for row in rows if row['covers_required'] == 'True')
    has_deletion_count = sum(1 for row in rows if row['full_deleted'] == 'True')
    
    # These counts should match the values in the YAML
    with open('test_output/test_output.yaml') as f:
        stats = yaml.safe_load(f)
    
    assert covers_required_count == stats['reads_covering_required']
    assert has_deletion_count == stats['reads_covered_full_deleted']

