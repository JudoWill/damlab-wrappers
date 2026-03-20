import os
from Bio import SeqIO


def test_bayesian_output_exists():
    assert os.path.exists("test_output/bayesian.consensus.fa")


def test_simple_output_exists():
    assert os.path.exists("test_output/simple.consensus.fa")


def test_bayesian_is_valid_fasta():
    records = list(SeqIO.parse("test_output/bayesian.consensus.fa", "fasta"))
    assert len(records) >= 1, "Bayesian consensus produced no sequences"


def test_simple_is_valid_fasta():
    records = list(SeqIO.parse("test_output/simple.consensus.fa", "fasta"))
    assert len(records) >= 1, "Simple consensus produced no sequences"


def test_bayesian_sequence_has_content():
    records = list(SeqIO.parse("test_output/bayesian.consensus.fa", "fasta"))
    assert len(records[0].seq) > 0, "Bayesian consensus sequence is empty"


def test_simple_sequence_has_content():
    records = list(SeqIO.parse("test_output/simple.consensus.fa", "fasta"))
    assert len(records[0].seq) > 0, "Simple consensus sequence is empty"


def test_max_reads_output_exists():
    assert os.path.exists("test_output/max_reads.consensus.fa")


def test_max_reads_is_valid_fasta():
    records = list(SeqIO.parse("test_output/max_reads.consensus.fa", "fasta"))
    assert len(records) >= 1, "max_reads consensus produced no sequences"


def test_max_reads_sequence_has_content():
    records = list(SeqIO.parse("test_output/max_reads.consensus.fa", "fasta"))
    assert len(records[0].seq) > 0, "max_reads consensus sequence is empty"
