import os


def test_compare_output_exists():
    """CRISPRessoCompare ran successfully and produced an output directory."""
    assert os.path.isdir("test_output/CRISPRessoCompare_sample1_vs_sample2")


def test_compare_log_exists():
    """CRISPRessoCompare running log exists."""
    assert os.path.exists(
        "test_output/CRISPRessoCompare_sample1_vs_sample2/CRISPRessoCompare_RUNNING_LOG.txt"
    )
