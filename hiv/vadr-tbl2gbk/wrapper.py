"""Wrapper for NCBI table2asn — convert a VADR feature table to GenBank format.

Takes the passing sequences FASTA and the NCBI 5-column feature table produced
by `v-annotate.pl` (files ending in `.vadr.pass.fa` and `.vadr.pass.tbl`) and
runs `table2asn` to produce an annotated GenBank flat file (`.gbf`) and an
ASN.1 submission file (`.sqn`).

Reference: https://www.ncbi.nlm.nih.gov/genbank/table2asn/

Typical upstream output from hiv/vadr-annotate:
    vadr/{sample}/ {sample}.vadr.pass.fa
    vadr/{sample}/ {sample}.vadr.pass.tbl
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import os
import shutil
import tempfile
from snakemake.shell import shell  # type: ignore

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
# Inputs — two accepted modes:
#
#   (a) Directory mode: input.vadr_dir points to a vadr-annotate output
#       directory; this wrapper locates *.vadr.pass.fa and *.vadr.pass.tbl
#       inside it automatically.  This is preferred because it avoids
#       referencing files inside a Snakemake directory() output, which
#       would cause a MissingInputException at DAG-build time.
#
#   (b) Explicit mode: input.sequences and input.feature_table are given
#       directly (useful outside the standard pipeline).
# ---------------------------------------------------------------------------

import glob as _glob

vadr_dir: str = snakemake.input.get("vadr_dir", "")
sequences: str = snakemake.input.get("sequences", "")
feature_table: str = snakemake.input.get("feature_table", "")

if vadr_dir:
    # Directory mode — find the pass files automatically
    if not os.path.isdir(vadr_dir):
        raise FileNotFoundError(
            f"input.vadr_dir '{vadr_dir}' does not exist or is not a directory."
        )
    pass_fa_matches = _glob.glob(os.path.join(vadr_dir, "*.vadr.pass.fa"))
    pass_tbl_matches = _glob.glob(os.path.join(vadr_dir, "*.vadr.pass.tbl"))
    if not pass_fa_matches:
        raise FileNotFoundError(
            f"No *.vadr.pass.fa file found in '{vadr_dir}'. "
            "v-annotate.pl may not have run successfully."
        )
    if not pass_tbl_matches:
        raise FileNotFoundError(
            f"No *.vadr.pass.tbl file found in '{vadr_dir}'. "
            "v-annotate.pl may not have run successfully."
        )
    sequences = pass_fa_matches[0]
    feature_table = pass_tbl_matches[0]
else:
    # Explicit mode
    if not sequences and len(snakemake.input) >= 1:
        sequences = str(snakemake.input[0])
    if not sequences:
        raise ValueError(
            "Provide either input.vadr_dir (vadr-annotate output directory) "
            "or input.sequences + input.feature_table explicitly."
        )
    if not feature_table and len(snakemake.input) >= 2:
        feature_table = str(snakemake.input[1])
    if not feature_table:
        raise ValueError(
            "input.feature_table must be provided when not using input.vadr_dir."
        )
    if not os.path.exists(sequences):
        raise FileNotFoundError(f"Input sequences file '{sequences}' does not exist.")
    if not os.path.exists(feature_table):
        raise FileNotFoundError(f"Input feature table '{feature_table}' does not exist.")

# Warn if the feature table is empty (no passing sequences)
if os.path.getsize(feature_table) == 0:
    import warnings
    warnings.warn(
        f"Feature table '{feature_table}' is empty — no passing sequences to convert.",
        stacklevel=1,
    )

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

output_gbf: str = str(snakemake.output.get("gbf", snakemake.output[0]))
output_sqn: str = str(snakemake.output.get("sqn", ""))

# ---------------------------------------------------------------------------
# Optional parameters
# ---------------------------------------------------------------------------

extra: str = snakemake.params.get("extra", "")
organism: str = snakemake.params.get("organism", "")
taxid: str = str(snakemake.params.get("taxid", ""))  # kept for API compat; not passed to table2asn (unrecognised modifier)

# Organism goes into table2asn's -j source-modifier string.
# Note: the bioconda table2asn build does not recognise [taxon=...] as a
# modifier; supply only [organism=...] here and let table2asn resolve the
# taxid from its built-in lookup.  -T is a boolean switch and does NOT
# accept a taxid value.
organism_arg = f'-j "[organism={organism}]"' if organism else ""

# ---------------------------------------------------------------------------
# Logging — table2asn writes a lot to stderr; capture it
# ---------------------------------------------------------------------------

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# ---------------------------------------------------------------------------
# Helpers: split multi-FASTA and multi-table into per-sequence dicts.
# table2asn -V b only accepts a single sequence per invocation; we must
# run it once per sequence and then concatenate the outputs.
# ---------------------------------------------------------------------------

def _parse_fasta(path):
    """Return OrderedDict of {seqid: [header_line, seq_lines...]}."""
    from collections import OrderedDict
    seqs = OrderedDict()
    cur_id, cur_lines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = cur_lines
                cur_id = line[1:].split()[0]
                cur_lines = [line]
            elif cur_id is not None:
                cur_lines.append(line)
    if cur_id is not None:
        seqs[cur_id] = cur_lines
    return seqs


def _parse_tbl(path):
    """Return dict of {seqid: [header_line, feature_lines...]}."""
    tbls = {}
    cur_id, cur_lines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Feature "):
                if cur_id is not None:
                    tbls[cur_id] = cur_lines
                cur_id = line.split()[1]
                cur_lines = [line]
            elif cur_id is not None:
                cur_lines.append(line)
    if cur_id is not None:
        tbls[cur_id] = cur_lines
    return tbls


# ---------------------------------------------------------------------------
# Execute: one table2asn call per sequence, then concatenate outputs.
# ---------------------------------------------------------------------------

fasta_seqs = _parse_fasta(sequences)
tbl_blocks = _parse_tbl(feature_table)

if not fasta_seqs:
    raise ValueError(f"No sequences found in '{sequences}'.")

gbf_parts = []
sqn_parts = []

with tempfile.TemporaryDirectory() as tmpdir:
    for seqid, fasta_lines in fasta_seqs.items():
        seq_fsa = os.path.join(tmpdir, f"{seqid}.fsa")
        seq_tbl = os.path.join(tmpdir, f"{seqid}.tbl")
        seq_gbf = os.path.join(tmpdir, f"{seqid}.gbf")
        seq_sqn = os.path.join(tmpdir, f"{seqid}.sqn")

        with open(seq_fsa, "w") as fh:
            fh.write("\n".join(fasta_lines) + "\n")

        if seqid in tbl_blocks:
            with open(seq_tbl, "w") as fh:
                fh.write("\n".join(tbl_blocks[seqid]) + "\n")
            tbl_arg = f"-f {seq_tbl}"
        else:
            tbl_arg = ""

        shell(
            "table2asn"
            " -i {seq_fsa}"
            " {tbl_arg}"
            " -V b"          # produce both .gbf (GenBank flat-file) and .sqn (ASN.1)
            " {organism_arg}"
            " {extra}"
            f" >> {snakemake.log[0]} 2>&1"
        )

        if os.path.exists(seq_gbf):
            gbf_parts.append(seq_gbf)
        if os.path.exists(seq_sqn):
            sqn_parts.append(seq_sqn)

    if not gbf_parts:
        raise FileNotFoundError(
            "table2asn did not produce any .gbf files. "
            "Check the log for errors. The feature table may be empty or malformed."
        )

    # Concatenate per-sequence outputs into the final files
    with open(output_gbf, "w") as out_fh:
        for part in gbf_parts:
            with open(part) as in_fh:
                out_fh.write(in_fh.read())

    if output_sqn:
        if sqn_parts:
            with open(output_sqn, "w") as out_fh:
                for part in sqn_parts:
                    with open(part) as in_fh:
                        out_fh.write(in_fh.read())
        else:
            raise FileNotFoundError(
                f"table2asn did not produce any .sqn files expected at '{output_sqn}'."
            )
