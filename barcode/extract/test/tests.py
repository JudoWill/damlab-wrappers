import pysam
import os
import yaml

def test_output_exists():
    """Test that output BAM file was created"""
    assert os.path.exists('test_output/test_output.bam')

def test_output_tags():
    """Test that output BAM has correct tags"""
    with pysam.AlignmentFile('test_output/test_output.bam', 'rb') as bam:
        for read in bam:
            if read.query_sequence:
                assert read.has_tag('CR')  # Has barcode tag
                assert read.has_tag('OX')  # Has UMI tag
                
                # Check tag lengths
                assert len(read.get_tag('CR')) == 10  # Barcode length
                assert len(read.get_tag('OX')) == 36  # Combined UMI length 

def test_metrics_file_exists():
    """Test that metrics file was created"""
    assert os.path.exists('test_output/test_metrics.yaml')

def test_metrics_file():
    """Test that metrics file is valid"""
    with open('test_output/test_metrics.yaml', 'r') as f:
        metrics = yaml.safe_load(f)
    assert metrics is not None
    
    expected_metrics = {
        'sample_name': 'sample_name',
        'total_reads': 100,
        'extraction_counts': {
            'barcode': 67,
            'left_umi': 44,
            'right_umi': 64
        },
        'pairwise_counts': {
            'all_three': 31,
            'barcode_left': 39,
            'barcode_right': 58,
            'left_right': 32
        }
    }
    assert metrics == expected_metrics
