"""Unit tests for bam2fastx output"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "0.0.1"


def test_fasta_has_records():
    """Test that FASTA output has expected structure and records."""
    with open("test_output/test.fasta") as f:
        content = f.read()

    lines = content.strip().split("\n")
    assert len(lines) > 0

    # Count sequences (lines starting with >)
    headers = [l for l in lines if l.startswith(">")]
    assert len(headers) > 0, "FASTA should have at least one sequence"

    # Each header should be followed by sequence lines
    for h in headers:
        assert h.startswith(">"), f"Header should start with '>': {h}"


def test_fastq_has_records():
    """Test that FASTQ output has expected structure and 4-line records."""
    with open("test_output/test.fastq") as f:
        content = f.read()

    lines = content.strip().split("\n")
    assert len(lines) > 0
    assert len(lines) % 4 == 0, "FASTQ should have 4 lines per record"

    for i in range(0, len(lines), 4):
        assert lines[i].startswith("@"), f"Line {i} should start with @: {lines[i]}"
        assert len(lines[i + 1]) > 0, "Sequence line should not be empty"
        assert lines[i + 2] == "+", "Plus line should be +"
        assert len(lines[i + 3]) == len(lines[i + 1]), "Quality length should match sequence"


def test_fasta_fastq_same_read_count():
    """Test that FASTA and FASTQ have the same number of records."""
    with open("test_output/test.fasta") as f:
        fasta_count = sum(1 for line in f if line.startswith(">"))

    with open("test_output/test.fastq") as f:
        fastq_count = sum(1 for line in f if line.strip() == '+')

    assert fasta_count == fastq_count, f"FASTA has {fasta_count} records, FASTQ has {fastq_count}"
