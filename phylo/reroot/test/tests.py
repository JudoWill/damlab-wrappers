import pytest # type: ignore
import os
import dendropy # type: ignore

def test_rerooted_tree_exists():
    """Test rerooted tree exists"""
    assert os.path.exists('test_output/test_rerooted.newick')

def test_rerooted_tree_is_valid():
    """Test rerooted tree is valid newick"""
    tree = dendropy.Tree.get(
        path='test_output/test_rerooted.newick',
        schema="newick"
    )
    assert tree is not None

def test_correct_root():
    """Test tree is rerooted correctly"""
    tree = dendropy.Tree.get(
        path='test_output/test_rerooted.newick',
        schema="newick"
    )
    
    for node in tree.leaf_node_iter():
        assert 'Seq1' != node.taxon.label, 'Root found as node taxon!'