import pytest # type: ignore
import os
import pysam # type: ignore

def test_single_output_exists():
    """Test that single-input BAM file was created"""
    assert os.path.exists('test_output/test.bam')

def test_directory_output_exists():
    """Test that directory-input BAM file was created"""
    assert os.path.exists('test_output/test_dir.bam')

def test_bam_is_valid():
    """Test that single-input output is a valid BAM file"""
    try:
        bam = pysam.AlignmentFile('test_output/test.bam', 'rb', check_sq=False)
        # Try to read a record to verify it's a valid BAM
        next(bam)
        bam.close()
    except Exception as e:
        pytest.fail(f"BAM file is not valid: {str(e)}")

def test_dir_bam_is_valid():
    """Test that directory-input output is a valid BAM file"""
    try:
        bam = pysam.AlignmentFile('test_output/test_dir.bam', 'rb', check_sq=False)
        # Try to read a record to verify it's a valid BAM
        next(bam)
        bam.close()
    except Exception as e:
        pytest.fail(f"Directory BAM file is not valid: {str(e)}")

def test_bam_has_reads():
    """Test that single-input BAM contains at least one read"""
    bam = pysam.AlignmentFile('test_output/test.bam', 'rb', check_sq=False)
    read_count = sum(1 for _ in bam)
    bam.close()
    assert read_count > 0, "BAM file contains no reads"

def test_dir_bam_has_reads():
    """Test that directory-input BAM contains at least one read"""
    bam = pysam.AlignmentFile('test_output/test_dir.bam', 'rb', check_sq=False)
    read_count = sum(1 for _ in bam)
    bam.close()
    assert read_count > 0, "Directory BAM file contains no reads" 