import os
import yaml

def test_output_exists():
    """Test that output files were created"""
    assert os.path.exists('test_output/test_gel.png')
    assert os.path.exists('test_output/test_metrics.yaml')

def test_yaml_output_format():
    """Test that output YAML has correct structure and values"""
    with open('test_output/test_metrics.yaml') as f:
        metrics = yaml.safe_load(f)
    
    # Check required keys
    required_keys = {'mode', 'lanes'}
    assert set(metrics.keys()) == required_keys
    
    # Check mode
    assert metrics['mode'] == 'read_length'
    
    # Check lanes structure
    assert 'test_sample' in metrics['lanes']
    lane_metrics = metrics['lanes']['test_sample']
    
    # Check lane metrics structure
    required_lane_keys = {
        'total_reads',
        'mean_size',
        'median_size',
        'min_size',
        'max_size'
    }
    assert set(lane_metrics.keys()) == required_lane_keys
    
    # Check that metrics are non-negative
    assert lane_metrics['total_reads'] >= 0
    assert lane_metrics['mean_size'] >= 0
    assert lane_metrics['median_size'] >= 0
    assert lane_metrics['min_size'] >= 0
    assert lane_metrics['max_size'] >= 0


def test_metrics_consistency():
    """Test that YAML metrics are consistent with expected values"""
    with open('test_output/test_metrics.yaml') as f:
        metrics = yaml.safe_load(f)
    
    lane_metrics = metrics['lanes']['test_sample']
    
    # Check that max_size is not exceeded
    assert lane_metrics['max_size'] <= 10000
    
    # Check that min_size is not negative
    assert lane_metrics['min_size'] >= 0
    
    # Check that mean and median are within reasonable bounds
    if lane_metrics['total_reads'] > 0:
        assert lane_metrics['mean_size'] >= lane_metrics['min_size']
        assert lane_metrics['mean_size'] <= lane_metrics['max_size']
        assert lane_metrics['median_size'] >= lane_metrics['min_size']
        assert lane_metrics['median_size'] <= lane_metrics['max_size'] 