import pytest # type: ignore
import os

def test_tree_output_exists():
    """Test tree exists"""
    assert os.path.exists('test_output/test.newick')

def test_tree_is_valid():
    """Test tree is valid"""
    with open('test_output/test.newick', 'r') as f:
        tree = f.read()
    assert tree.strip() != ''