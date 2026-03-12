"""Tests for the hiv/vadr-build wrapper.

Checks that v-build.pl produces the expected output directory and key files
for accession K03455 (HIV-1 HXB2).
"""

import os
import pytest

MODEL_DIR = "test_output/K03455"
ACCESSION = "K03455"


def test_output_directory_exists():
    """The model directory must be created by v-build.pl."""
    assert os.path.isdir(MODEL_DIR), f"Model directory '{MODEL_DIR}' was not created."


def test_cm_file_exists():
    """A covariance model (.cm) file must be present."""
    cm = os.path.join(MODEL_DIR, f"{ACCESSION}.cm")
    assert os.path.exists(cm), f"Expected CM file '{cm}' not found."


def test_minfo_file_exists():
    """A model-info (.minfo) file must be present."""
    minfo = os.path.join(MODEL_DIR, f"{ACCESSION}.minfo")
    assert os.path.exists(minfo), f"Expected minfo file '{minfo}' not found."


def test_reference_fasta_exists():
    """The reference FASTA used to build the model must be present."""
    fa = os.path.join(MODEL_DIR, f"{ACCESSION}.fa")
    assert os.path.exists(fa), f"Expected reference FASTA '{fa}' not found."


def test_reference_fasta_nonempty():
    """The reference FASTA must contain sequence data."""
    fa = os.path.join(MODEL_DIR, f"{ACCESSION}.fa")
    if os.path.exists(fa):
        assert os.path.getsize(fa) > 0, f"Reference FASTA '{fa}' is empty."


def test_cm_file_nonempty():
    """The CM file must not be empty."""
    cm = os.path.join(MODEL_DIR, f"{ACCESSION}.cm")
    if os.path.exists(cm):
        assert os.path.getsize(cm) > 0, f"CM file '{cm}' is empty."


def test_blast_nucleotide_db_exists():
    """At least one BLAST nucleotide index file (.nhr/.nin/.nsq) must be present."""
    blast_files = [
        f for f in os.listdir(MODEL_DIR)
        if f.endswith((".nhr", ".nin", ".nsq"))
    ]
    assert len(blast_files) > 0, (
        f"No BLAST nucleotide DB files found in '{MODEL_DIR}'. "
        "Expected at least one .nhr/.nin/.nsq file."
    )
