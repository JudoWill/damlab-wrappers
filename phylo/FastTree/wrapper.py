"""Wrapper for FastTree phylogenetic tree inference.

This wrapper provides a Snakemake interface to FastTree, a program for inferring
approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide
or protein sequences.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.1"

import os
from typing import Optional

from snakemake.shell import shell  # type: ignore

# This is a common pattern in Snakemake wrappers
# It allows the wrapper to be imported without snakemake being in the global namespace
# This is useful for testing and linting
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Check if version is specified and compatible
if hasattr(snakemake, "params") and "version" in snakemake.params:
    requested_version = snakemake.params.version
    if requested_version != __version__:
        print(f"Warning: Requested version {requested_version} does not match wrapper version {__version__}")

# Get input and output files
alignment: str = snakemake.input[0]
output: str = snakemake.output[0]

# Validate input file exists
if not os.path.exists(alignment):
    raise FileNotFoundError(f"Input file {alignment} does not exist")

# Get parameters with defaults
# This is a common pattern in wrappers - providing sensible defaults
extra: str = snakemake.params.get("extra", "")
gtr: bool = snakemake.params.get("gtr", False)
gamma: bool = snakemake.params.get("gamma", False)
boot: Optional[int] = snakemake.params.get("boot", None)
nt: bool = snakemake.params.get("nt", True)  # Default to nucleotide sequences

# Build command line arguments
cmd_args = []

# Add model parameters
if gtr:
    cmd_args.append("-gtr")
if gamma:
    cmd_args.append("-gamma")
if boot is not None:
    cmd_args.extend(["-boot", str(boot)])
if nt:
    cmd_args.append("-nt")

# Add input and output
cmd_args.extend([alignment])

# Add extra parameters if specified
if extra:
    cmd_args.append(extra)

# Setup logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Execute FastTree
shell(f"FastTree {' '.join(cmd_args)} > {output} {log}")