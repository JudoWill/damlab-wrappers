"""Wrapper for VADR v-build.pl — build a homology model from a GenBank accession.

`v-build.pl` fetches a RefSeq/GenBank accession and creates VADR model files
(CM, BLAST databases, feature info) used by `v-annotate.pl` to validate and
annotate similar sequences.  Only sequences ≤ 25 kb are supported due to
memory requirements.

By default the accession is fetched from NCBI.  Two optional inputs allow
fully offline operation:

  input.fasta   — local FASTA file (passed to v-build.pl with --infa)
  input.genbank — local GenBank file (passed with --gb --ingb); this replaces
                  both the feature table fetch and the FASTA fetch when
                  combined with input.fasta.

When using local files downloaded from NCBI (e.g. K03455.gb / K03455.fasta),
run them through the `hiv/vadr-genbank` wrapper first to strip version suffixes
(.1, .2 …) that cause VADR's accession-matching checks to fail.

Reference: https://github.com/ncbi/vadr/blob/master/documentation/build.md
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
# Required parameters
# ---------------------------------------------------------------------------

accession: str = snakemake.params.get("accession", "")
if not accession:
    raise ValueError("params.accession must be provided (e.g. 'K03455')")

# ---------------------------------------------------------------------------
# Output — must be a single directory
# ---------------------------------------------------------------------------

if len(snakemake.output) != 1:
    raise ValueError("Exactly one output must be specified (the model directory to create).")

output_dir: str = str(snakemake.output[0])

# Snakemake pre-creates output directories; v-build.pl expects to create the
# directory itself, so we remove it if it was pre-created and is still empty.
if os.path.isdir(output_dir) and not os.listdir(output_dir):
    os.rmdir(output_dir)
elif os.path.exists(output_dir):
    raise FileExistsError(
        f"Output directory '{output_dir}' already exists and is non-empty. "
        "Remove it before re-running."
    )

# ---------------------------------------------------------------------------
# Optional local file inputs (offline mode)
# ---------------------------------------------------------------------------

local_fasta: str = snakemake.input.get("fasta", "")
local_genbank: str = snakemake.input.get("genbank", "")

if local_fasta and not os.path.exists(local_fasta):
    raise FileNotFoundError(f"input.fasta '{local_fasta}' does not exist.")
if local_genbank and not os.path.exists(local_genbank):
    raise FileNotFoundError(f"input.genbank '{local_genbank}' does not exist.")

infa_arg = f"--infa {local_fasta}" if local_fasta else ""
ingb_arg = f"--gb --ingb {local_genbank}" if local_genbank else ""

# ---------------------------------------------------------------------------
# Optional parameters
# ---------------------------------------------------------------------------

extra: str = snakemake.params.get("extra", "")
group: str = snakemake.params.get("group", "")
subgroup: str = snakemake.params.get("subgroup", "")

group_arg = f"--group {group}" if group else ""
subgroup_arg = f"--subgroup {subgroup}" if subgroup else ""

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# ---------------------------------------------------------------------------
# Execute
# ---------------------------------------------------------------------------

shell(
    "v-build.pl"
    " {infa_arg}"
    " {ingb_arg}"
    " {group_arg}"
    " {subgroup_arg}"
    " {extra}"
    " {accession}"
    " {output_dir}"
    " {log}"
)
