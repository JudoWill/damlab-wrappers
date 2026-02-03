import os

def test_alignment_output_exists():
    """Test alignment exists"""
    assert os.path.exists('test_output/test.fasta')

def test_alignment_is_valid():
    """Test alignment is valid"""
    with open('test_output/test.fasta', 'r') as f:
        alignment = f.read()
    assert alignment.strip() != '' 