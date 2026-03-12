"""Integration tests for the VADR annotation pipeline.

Verifies that vadr-build → vadr-annotate → vadr-tbl2gbk produces
well-formed, non-empty outputs for each accession in ACCESSIONS, using
the bundled haplotypes.fa test sequences.
"""

import os
import pytest

# Must match the ACCESSIONS list in the Snakefile
ACCESSIONS = [
    "K03455",
    "AF324493",
]


# ---------------------------------------------------------------------------
# Parametrised fixtures — one test instance per accession
# ---------------------------------------------------------------------------

@pytest.fixture(params=ACCESSIONS)
def accession(request):
    return request.param


def model_dir(accession):
    return f"test_output/model/{accession}"

def annotation_dir(accession):
    return f"test_output/{accession}/annotation"

def gbf_path(accession):
    return f"test_output/{accession}/haplotypes.gbf"

def sqn_path(accession):
    return f"test_output/{accession}/haplotypes.sqn"


# ---------------------------------------------------------------------------
# vadr-build outputs
# ---------------------------------------------------------------------------

def test_model_directory_exists(accession):
    mdir = model_dir(accession)
    assert os.path.isdir(mdir), f"[{accession}] Model directory '{mdir}' was not created."

def test_model_cm_exists(accession):
    cm = os.path.join(model_dir(accession), f"{accession}.vadr.cm")
    assert os.path.exists(cm), f"[{accession}] CM file not found: {cm}"

def test_model_minfo_exists(accession):
    minfo = os.path.join(model_dir(accession), f"{accession}.vadr.minfo")
    assert os.path.exists(minfo), f"[{accession}] minfo file not found: {minfo}"

def test_model_blast_db_exists(accession):
    mdir = model_dir(accession)
    nhr_files = [f for f in os.listdir(mdir) if f.endswith(".nhr")] if os.path.isdir(mdir) else []
    assert nhr_files, f"[{accession}] No BLAST nucleotide .nhr file found in '{mdir}'."


# ---------------------------------------------------------------------------
# vadr-annotate outputs
# ---------------------------------------------------------------------------

def _annotation_files(accession, suffix):
    adir = annotation_dir(accession)
    if not os.path.isdir(adir):
        return []
    return [f for f in os.listdir(adir) if f.endswith(suffix)]

def test_annotation_directory_exists(accession):
    adir = annotation_dir(accession)
    assert os.path.isdir(adir), f"[{accession}] Annotation directory '{adir}' was not created."

def test_annotation_sqc_exists(accession):
    assert _annotation_files(accession, ".sqc"), f"[{accession}] No .sqc classification file found."

def test_annotation_ftr_exists(accession):
    assert _annotation_files(accession, ".ftr"), f"[{accession}] No .ftr feature file found."

def test_annotation_pass_fa_exists(accession):
    assert _annotation_files(accession, ".pass.fa"), f"[{accession}] No .pass.fa — no sequences passed."

def test_annotation_pass_tbl_nonempty(accession):
    files = _annotation_files(accession, ".pass.tbl")
    if not files:
        pytest.skip(f"[{accession}] No .pass.tbl found.")
    path = os.path.join(annotation_dir(accession), files[0])
    assert os.path.getsize(path) > 0, f"[{accession}] .pass.tbl '{path}' is empty."

def test_annotation_has_passing_sequences(accession):
    """At least one haplotype must pass VADR annotation."""
    files = _annotation_files(accession, ".pass.fa")
    if not files:
        pytest.skip(f"[{accession}] No .pass.fa found.")
    path = os.path.join(annotation_dir(accession), files[0])
    assert os.path.getsize(path) > 0, (
        f"[{accession}] .pass.fa is empty — all haplotypes failed VADR."
    )


# ---------------------------------------------------------------------------
# vadr-tbl2gbk outputs
# ---------------------------------------------------------------------------

def test_gbf_exists(accession):
    gbf = gbf_path(accession)
    assert os.path.exists(gbf), f"[{accession}] GenBank file '{gbf}' was not created."

def test_sqn_exists(accession):
    sqn = sqn_path(accession)
    assert os.path.exists(sqn), f"[{accession}] ASN.1 file '{sqn}' was not created."

def test_gbf_nonempty(accession):
    gbf = gbf_path(accession)
    assert os.path.getsize(gbf) > 0, f"[{accession}] GenBank file '{gbf}' is empty."

def test_gbf_has_locus(accession):
    with open(gbf_path(accession)) as fh:
        content = fh.read()
    assert "LOCUS" in content, f"[{accession}] GenBank file missing LOCUS header."

def test_gbf_has_features(accession):
    with open(gbf_path(accession)) as fh:
        content = fh.read()
    assert "FEATURES" in content, f"[{accession}] GenBank file missing FEATURES section."

def test_gbf_has_origin(accession):
    with open(gbf_path(accession)) as fh:
        content = fh.read()
    assert "ORIGIN" in content, f"[{accession}] GenBank file missing ORIGIN section."

def test_gbf_properly_terminated(accession):
    with open(gbf_path(accession)) as fh:
        content = fh.read()
    assert "//" in content, f"[{accession}] GenBank file missing '//' record terminator."

def test_gbf_organism_annotation(accession):
    """Source feature should carry the organism name passed to table2asn."""
    with open(gbf_path(accession)) as fh:
        content = fh.read()
    assert "immunodeficiency" in content.lower(), (
        f"[{accession}] Expected organism annotation not found in GenBank file."
    )
