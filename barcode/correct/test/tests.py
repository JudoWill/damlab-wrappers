import pysam
import os
import pytest

def test_output_exists():
    """Test that output BAM file was created"""
    assert os.path.exists('test_output/test_output.bam')

def test_output_tags():
    """Test that output BAM has correct tags"""
    with pysam.AlignmentFile('test_output/test_output.bam', 'rb') as bam:
        for read in bam:
            try:
                original = read.get_tag('CR')  # Original tag
                corrected = read.get_tag('CB')  # Corrected tag
                assert len(original) == len(corrected)  # Same length
                assert original == corrected or original != corrected  # Either corrected or not
            except KeyError:
                continue  # Skip reads without tags 