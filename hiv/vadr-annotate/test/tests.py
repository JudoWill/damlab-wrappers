"""Tests for the hiv/vadr-annotate wrapper.

Checks that v-annotate.pl produces the expected output directory and key
output files for the bundled K03455 test sequence.
"""

import os
import pytest

ANNOTATION_DIR = "test_output/annotation"


def test_output_directory_exists():
    """The annotation directory must be created by v-annotate.pl."""
    assert os.path.isdir(ANNOTATION_DIR), (
        f"Annotation directory '{ANNOTATION_DIR}' was not created."
    )


def _find_files_with_suffix(suffix):
    """Return all files in ANNOTATION_DIR with the given suffix."""
    if not os.path.isdir(ANNOTATION_DIR):
        return []
    return [
        f for f in os.listdir(ANNOTATION_DIR)
        if f.endswith(suffix)
    ]


def test_sqc_file_exists():
    """A per-sequence classification summary (.sqc) file must be present."""
    files = _find_files_with_suffix(".sqc")
    assert len(files) > 0, (
        f"No .sqc file found in '{ANNOTATION_DIR}'. "
        "v-annotate.pl may have failed."
    )


def test_ftr_file_exists():
    """A per-feature annotation table (.ftr) file must be present."""
    files = _find_files_with_suffix(".ftr")
    assert len(files) > 0, (
        f"No .ftr file found in '{ANNOTATION_DIR}'."
    )


def test_alt_file_exists():
    """An alert detail listing (.alt) file must be present."""
    files = _find_files_with_suffix(".alt")
    assert len(files) > 0, (
        f"No .alt file found in '{ANNOTATION_DIR}'."
    )


def test_log_file_exists():
    """A run log (.log) file must be present in the annotation directory."""
    files = _find_files_with_suffix(".log")
    assert len(files) > 0, (
        f"No .log file found in '{ANNOTATION_DIR}'."
    )


def test_output_fasta_exists():
    """At least one output FASTA (.pass.fa or .fail.fa) must be present."""
    pass_files = _find_files_with_suffix(".pass.fa")
    fail_files = _find_files_with_suffix(".fail.fa")
    assert len(pass_files) + len(fail_files) > 0, (
        f"No .pass.fa or .fail.fa files found in '{ANNOTATION_DIR}'. "
        "v-annotate.pl did not produce sequence output."
    )


def test_sqc_file_nonempty():
    """The .sqc classification summary must not be empty."""
    files = _find_files_with_suffix(".sqc")
    if files:
        path = os.path.join(ANNOTATION_DIR, files[0])
        assert os.path.getsize(path) > 0, f".sqc file '{path}' is empty."


def test_sqc_contains_sequence_entry():
    """The .sqc file must contain at least one sequence classification row."""
    files = _find_files_with_suffix(".sqc")
    if not files:
        pytest.skip("No .sqc file found; skipping content check.")
    path = os.path.join(ANNOTATION_DIR, files[0])
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#") and l.strip()]
    assert len(lines) > 0, (
        f".sqc file '{path}' contains no non-comment lines."
    )
