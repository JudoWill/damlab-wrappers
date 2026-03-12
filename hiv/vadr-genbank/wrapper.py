"""Wrapper: normalise GenBank / FASTA files for use with VADR v-build.pl.

VADR's GenBank parser compares the VERSION line in a `.gb` file against the
bare accession argument passed to `v-build.pl`.  Files downloaded directly
from NCBI contain versioned accessions (e.g. ``VERSION K03455.1``,
``>K03455.1``), which cause v-build.pl to abort with:

    ERROR in sqf_GenbankParse … version/accession mismatch for K03455,
    line: VERSION K03455.1

This wrapper rewrites the offending fields so the bare accession (no ``.N``
suffix) appears wherever VADR looks for it:

* GenBank flat file — ``LOCUS``, ``ACCESSION``, and ``VERSION`` lines
* FASTA file — the accession token (first whitespace-delimited word) in every
  ``>`` header line

The original files are **never modified**.  Output files are written to the
paths specified in ``output.genbank`` / ``output.fasta``.

If a file does not contain any versioned form of the accession the content is
copied unchanged, so this wrapper is safe to use on files that are already
clean or on custom references with no version suffix.

Reference: https://github.com/ncbi/vadr/blob/master/documentation/build.md
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import os
import re
import shutil


if "snakemake" not in locals():
    import snakemake  # type: ignore

if hasattr(snakemake, "params") and "version" in snakemake.params:
    requested_version = snakemake.params.version
    if requested_version != __version__:
        print(
            f"Warning: Requested version {requested_version} does not match "
            f"wrapper version {__version__}"
        )

# ---------------------------------------------------------------------------
# Required parameter
# ---------------------------------------------------------------------------

accession: str = snakemake.params.get("accession", "")
if not accession:
    raise ValueError(
        "params.accession must be provided (e.g. 'K03455').  "
        "It is the bare accession — no version suffix."
    )

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

input_genbank: str = snakemake.input.get("genbank", "")
input_fasta: str = snakemake.input.get("fasta", "")

if not input_genbank and not input_fasta:
    raise ValueError("At least one of input.genbank or input.fasta must be provided.")

if input_genbank and not os.path.exists(input_genbank):
    raise FileNotFoundError(f"input.genbank '{input_genbank}' does not exist.")
if input_fasta and not os.path.exists(input_fasta):
    raise FileNotFoundError(f"input.fasta '{input_fasta}' does not exist.")

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

output_genbank: str = snakemake.output.get("genbank", "")
output_fasta: str = snakemake.output.get("fasta", "")

if input_genbank and not output_genbank:
    raise ValueError("output.genbank must be provided when input.genbank is given.")
if input_fasta and not output_fasta:
    raise ValueError("output.fasta must be provided when input.fasta is given.")

# ---------------------------------------------------------------------------
# Helper — ensure parent directory exists
# ---------------------------------------------------------------------------

def _ensure_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


# ---------------------------------------------------------------------------
# GenBank normalisation
#
# VADR's sqf_GenbankParse (in sequip/sqp_seqfile.pm) uses the first token of
# the LOCUS line as the primary accession key ($acc), then validates that the
# VERSION line (after stripping its own ".N" suffix) matches $acc exactly.
#
# Two problems can arise with files downloaded from NCBI:
#
#   1. Legacy LOCUS name differs from the accession.
#      e.g. K03455.gb has  LOCUS HIVHXB2CG  but VERSION K03455.1
#      VADR sets $acc = "HIVHXB2CG", strips VERSION to "K03455", and
#      "K03455" ne "HIVHXB2CG" → "version/accession mismatch" error.
#
#   2. VERSION has a ".N" version suffix.
#      Even when LOCUS == accession, VADR's strip+compare expects the bare
#      accession.  The strip is applied by seq_StripVersion, so "K03455.1"
#      → "K03455" and the comparison succeeds — but only after the LOCUS
#      check above also passes.
#
# Fix: rewrite the LOCUS name to the bare accession AND strip the ".N"
# version suffix from the VERSION line.
# ---------------------------------------------------------------------------

if input_genbank:
    with open(input_genbank) as fh:
        gb_content = fh.read()

    changes = []
    normalised = gb_content

    # 1. Replace the LOCUS name token with the bare accession.
    #    LOCUS line format: "LOCUS       <name>   <len> bp ..."
    #    We only replace when the name is NOT already the bare accession.
    _locus_re = re.compile(
        r"^(LOCUS\s+)(\S+)(\s+\d+\s+)",
        re.MULTILINE,
    )
    def _fix_locus(m):
        old_name = m.group(2)
        if old_name != accession:
            changes.append(f"LOCUS name '{old_name}' → '{accession}'")
        return m.group(1) + accession + m.group(3)

    normalised = _locus_re.sub(_fix_locus, normalised)

    # 2. Strip ".N" version suffix from the VERSION line.
    _ver_re = re.compile(
        r"^(VERSION\s+)" + re.escape(accession) + r"(\.\d+)",
        re.MULTILINE,
    )
    n_ver = len(_ver_re.findall(normalised))
    if n_ver:
        normalised = _ver_re.sub(lambda m: m.group(1) + accession, normalised)
        changes.append(f"VERSION suffix stripped ({n_ver} occurrence(s))")

    _ensure_dir(output_genbank)
    with open(output_genbank, "w") as fh:
        fh.write(normalised)

    if changes:
        print(
            f"[vadr-genbank] GenBank: normalised {os.path.basename(input_genbank)}: "
            + "; ".join(changes)
        )
    else:
        print(
            f"[vadr-genbank] GenBank: no changes needed in "
            f"{os.path.basename(input_genbank)} — copied unchanged"
        )

# ---------------------------------------------------------------------------
# FASTA — copy unchanged
#
# v-build.pl requires FASTA headers to be in "<accession>.<N>" format
# (e.g. ">K03455.1").  The Perl regex is:
#
#   if ($mdl_name_ver =~ /(\S+)\.\d+/) { check $1 eq $mdl_name }
#   else { ERROR "expected accession.version …" }
#
# Stripping the ".N" suffix causes an unconditional failure in that branch,
# so the FASTA must be passed to v-build.pl with the version suffix intact.
# We simply copy it to the output path.
# ---------------------------------------------------------------------------

if input_fasta:
    _ensure_dir(output_fasta)
    shutil.copy2(input_fasta, output_fasta)
    print(
        f"[vadr-genbank] FASTA: copied unchanged "
        f"(version suffix must be kept for v-build.pl)"
    )
