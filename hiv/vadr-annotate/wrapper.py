"""Wrapper for VADR v-annotate.pl — classify and annotate viral sequences.

`v-annotate.pl` uses VADR model files (created by `v-build.pl`) to validate
and annotate sequences in a FASTA file.  It reports passing/failing sequences
and generates feature tables, alignment files, and detailed alert listings.

Reference: https://github.com/ncbi/vadr/blob/master/documentation/annotate.md
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import os
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
# Required input: sequences FASTA
# ---------------------------------------------------------------------------

sequences: str = snakemake.input.get("sequences", "")
if not sequences:
    # Fall back to positional input[0]
    if len(snakemake.input) >= 1:
        sequences = str(snakemake.input[0])
    else:
        raise ValueError(
            "An input sequences FASTA file must be provided via input.sequences "
            "or as the first positional input."
        )

if not os.path.exists(sequences):
    raise FileNotFoundError(f"Input sequences file '{sequences}' does not exist.")

# ---------------------------------------------------------------------------
# Optional input: model library directory (from vadr-build)
# ---------------------------------------------------------------------------

mdir: str = snakemake.input.get("mdir", "")
if mdir and not os.path.isdir(mdir):
    raise FileNotFoundError(
        f"Model library directory '{mdir}' does not exist. "
        "Run hiv/vadr-build first, or point to an existing VADR model library."
    )

# ---------------------------------------------------------------------------
# Output — must be a single directory
# ---------------------------------------------------------------------------

if len(snakemake.output) != 1:
    raise ValueError(
        "Exactly one output must be specified (the annotation directory to create)."
    )

output_dir: str = str(snakemake.output[0])

# v-annotate.pl expects to create the directory itself.
if os.path.isdir(output_dir) and not os.listdir(output_dir):
    os.rmdir(output_dir)
elif os.path.exists(output_dir):
    raise FileExistsError(
        f"Output directory '{output_dir}' already exists and is non-empty. "
        "Remove it before re-running."
    )

# ---------------------------------------------------------------------------
# Optional parameters
# ---------------------------------------------------------------------------

extra: str = snakemake.params.get("extra", "")
noseqnamemax: bool = snakemake.params.get("noseqnamemax", False)

# All VADR alert codes that are *fatal by default* and can be downgraded to
# non-fatal via --alt_pass.  The only alert that cannot be changed is
# `ftskipfl` (always-fatal).
# Source: VADR annotate.md §"Description of alerts that are fatal by default"
# https://github.com/NLM-DIR/vadr/blob/master/documentation/annotate.md
_HIV_ALT_PASS = (
    # sequence-level
    "incsbgrp,incgroup,lowcovrg,dupregin,discontn,indfstrn,"
    "lowsim5s,lowsim3s,lowsimis,nmiscftr,deletins,"
    # feature-level
    "mutstart,mutendcd,mutendns,mutendex,unexleng,"
    "cdsstopn,cdsstopp,fsthicft,fsthicfi,fstukcft,fstukcfi,"
    "mutspst5,mutspst3,peptrans,pepadjcy,"
    "indfantp,indfantn,"
    "indf5gap,indf5lcn,indf5plg,indf5pst,"
    "indf3gap,indf3lcn,indf3plg,indf3pst,"
    "indfstrp,insertnp,deletinp,deletinf,"
    "lowsim5n,lowsim5l,lowsim3n,lowsim3l,lowsimin,lowsimil"
)

mode: str = snakemake.params.get("mode", "")
if mode == "hiv":
    alt_pass_arg = f"--alt_pass {_HIV_ALT_PASS}"
elif mode and mode != "hiv":
    raise ValueError(f"Unknown mode '{mode}'. Supported modes: 'hiv'.")
else:
    alt_pass_arg = ""

# --mdir and --mkey handling.
# When --mdir points to a single-model directory produced by vadr-build, the
# .minfo file is named "{accession}.vadr.minfo" (e.g. "K03455.vadr.minfo").
# v-annotate.pl's default --mkey ("calici") won't find it, so we detect the
# key automatically from the only .minfo file present in the directory.
mdir_arg = ""
mkey_arg = ""
if mdir:
    import glob as _glob
    minfo_files = _glob.glob(os.path.join(mdir, "*.minfo"))
    if len(minfo_files) == 1:
        mkey = os.path.basename(minfo_files[0]).removesuffix(".minfo")
        mkey_arg = f"--mkey {mkey}"
    elif len(minfo_files) > 1:
        # Multiple models — let the user specify via extra or rely on the default
        import warnings
        warnings.warn(
            f"Multiple .minfo files found in '{mdir}'; not setting --mkey automatically. "
            "Use params.extra to pass '--mkey <key>' if needed.",
            stacklevel=1,
        )
    mdir_arg = f"--mdir {mdir}"

threads_arg = (
    f"--split --cpu {snakemake.threads}" if snakemake.threads > 1 else ""
)
noseqnamemax_arg = "--noseqnamemax" if noseqnamemax else ""

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# ---------------------------------------------------------------------------
# Execute
# ---------------------------------------------------------------------------

shell(
    "v-annotate.pl"
    " {mdir_arg}"
    " {mkey_arg}"
    " {alt_pass_arg}"
    " {threads_arg}"
    " {noseqnamemax_arg}"
    " {extra}"
    " {sequences}"
    " {output_dir}"
    " {log}"
)
