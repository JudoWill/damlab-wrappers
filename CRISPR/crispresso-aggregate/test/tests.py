import os


def test_aggregate_output_exists():
    """CRISPRessoAggregate ran successfully and produced an output directory."""
    assert os.path.isdir("test_output/CRISPRessoAggregate_on_all_samples")


def test_aggregate_log_exists():
    """CRISPRessoAggregate running log exists."""
    assert os.path.exists(
        "test_output/CRISPRessoAggregate_on_all_samples/CRISPRessoAggregate_RUNNING_LOG.txt"
    )
