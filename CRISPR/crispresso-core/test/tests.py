import os


def test_param_amplicon_output_exists():
    """CRISPResso ran successfully with amplicon provided as a param string."""
    assert os.path.isdir("test_output/CRISPResso_on_param_amplicon")


def test_param_amplicon_log_exists():
    """CRISPResso running log exists for param-amplicon run."""
    assert os.path.exists(
        "test_output/CRISPResso_on_param_amplicon/CRISPResso_RUNNING_LOG.txt"
    )


def test_param_amplicon_quantification_exists():
    """CRISPResso quantification output exists for param-amplicon run."""
    assert os.path.exists(
        "test_output/CRISPResso_on_param_amplicon/"
        "CRISPResso_quantification_of_editing_frequency.txt"
    )


def test_fasta_amplicon_output_exists():
    """CRISPResso ran successfully with amplicon provided as a FASTA file."""
    assert os.path.isdir("test_output/CRISPResso_on_fasta_amplicon")


def test_fasta_amplicon_log_exists():
    """CRISPResso running log exists for fasta-amplicon run."""
    assert os.path.exists(
        "test_output/CRISPResso_on_fasta_amplicon/CRISPResso_RUNNING_LOG.txt"
    )


def test_fasta_amplicon_quantification_exists():
    """CRISPResso quantification output exists for fasta-amplicon run."""
    assert os.path.exists(
        "test_output/CRISPResso_on_fasta_amplicon/"
        "CRISPResso_quantification_of_editing_frequency.txt"
    )
