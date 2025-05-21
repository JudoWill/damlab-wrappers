import pysam
import yaml

def test_inner_output_tags():
    """Test that output BAM has correct tags"""
    with pysam.AlignmentFile('test_output/test_inner.bam', 'rb') as bam:
        for read in bam:
            if read.query_sequence:
                assert read.has_tag('CR')  # Has barcode tag
                assert read.has_tag('OX')  # Has UMI tag
                
                # Check tag lengths
                assert len(read.get_tag('CR')) in {0, 10}  # Barcode length
                assert len(read.get_tag('OX')) in {0, 18, 36}  # Combined UMI length 

def test_whole_output_tags():
    """Test that output BAM has correct tags"""
    with pysam.AlignmentFile('test_output/test_whole.bam', 'rb') as bam:
        for read in bam:
            if read.query_sequence:
                assert read.has_tag('CR')  # Has barcode tag
                assert read.has_tag('OX')  # Has UMI tag

                # Check tag lengths
                assert len(read.get_tag('CR')) in {0, 34}  # Barcode length
                assert len(read.get_tag('OX')) in {0, 18, 36}  # Combined UMI length 



def test_inner_metrics_file():
    """Test that inner metrics file is valid"""
    with open('test_output/test_inner.yaml', 'r') as f:
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

def test_whole_metrics_file():
    """Test that whole metrics file is valid"""
    with open('test_output/test_whole.yaml', 'r') as f:
        metrics = yaml.safe_load(f)
    assert metrics is not None
    
    expected_metrics = {
        'sample_name': 'sample_name',
        'total_reads': 100,
        'extraction_counts': {
            'barcode': 69,
            'left_umi': 44,
            'right_umi': 64
        },
        'pairwise_counts': {
            'all_three': 32,
            'barcode_left': 40,
            'barcode_right': 60,
            'left_right': 32
        }
    }
    assert metrics == expected_metrics