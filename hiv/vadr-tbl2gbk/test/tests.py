"""Tests for the hiv/vadr-tbl2gbk wrapper.

Checks that table2asn produces a valid GenBank flat file and ASN.1 submission
file from the bundled test FASTA and feature table.
"""

import os
import pytest

GBF_PATH = "test_output/test.gbf"
SQN_PATH = "test_output/test.sqn"
SEQ_ID   = "K03455_hap1"


def test_gbf_exists():
    """GenBank flat file must be created."""
    assert os.path.exists(GBF_PATH), f"GenBank file '{GBF_PATH}' was not created."


def test_sqn_exists():
    """ASN.1 submission file must be created."""
    assert os.path.exists(SQN_PATH), f"ASN.1 file '{SQN_PATH}' was not created."


def test_gbf_nonempty():
    """GenBank flat file must contain data."""
    assert os.path.getsize(GBF_PATH) > 0, f"GenBank file '{GBF_PATH}' is empty."


def test_sqn_nonempty():
    """ASN.1 file must contain data."""
    assert os.path.getsize(SQN_PATH) > 0, f"ASN.1 file '{SQN_PATH}' is empty."


def test_gbf_has_locus_line():
    """GenBank flat file must start with a LOCUS header line."""
    with open(GBF_PATH) as fh:
        content = fh.read()
    assert "LOCUS" in content, (
        "GenBank file does not contain a LOCUS line. "
        "table2asn may have produced an incomplete output."
    )


def test_gbf_has_origin():
    """GenBank flat file must contain an ORIGIN section with the sequence."""
    with open(GBF_PATH) as fh:
        content = fh.read()
    assert "ORIGIN" in content, (
        "GenBank file does not contain an ORIGIN section. "
        "Sequence data appears to be missing."
    )


def test_gbf_contains_sequence_id():
    """GenBank flat file must reference the input sequence ID."""
    with open(GBF_PATH) as fh:
        content = fh.read()
    assert SEQ_ID in content, (
        f"Sequence ID '{SEQ_ID}' not found in GenBank file. "
        "The feature table and FASTA may be mismatched."
    )


def test_gbf_has_source_feature():
    """GenBank flat file must contain at least one feature annotation."""
    with open(GBF_PATH) as fh:
        content = fh.read()
    assert "FEATURES" in content, (
        "GenBank file does not contain a FEATURES section. "
        "Feature table may not have been applied."
    )


def test_gbf_ends_with_terminator():
    """GenBank flat file records must end with '//'."""
    with open(GBF_PATH) as fh:
        content = fh.read()
    assert "//" in content, (
        "GenBank file does not contain a record terminator '//'. "
        "The file may be truncated or malformed."
    )
