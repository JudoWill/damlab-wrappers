"""Wrapper for MUSCLE multiple sequence alignment tool.

This wrapper provides a Snakemake interface to the MUSCLE multiple sequence alignment tool.
It supports all standard MUSCLE parameters through the 'extra' parameter and provides
sensible defaults for common use cases.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.1"

import os
from pathlib import Path
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
input_seqs: str = snakemake.input[0]
output_aln: str = snakemake.output[0]

# Validate input file exists
if not os.path.exists(input_seqs):
    raise FileNotFoundError(f"Input file {input_seqs} does not exist")

# Get parameters with defaults
# This is a common pattern in wrappers - providing sensible defaults
extra: str = snakemake.params.get("extra", "")
threads: int = snakemake.threads
maxiters: Optional[int] = snakemake.params.get("maxiters", None)
diags: bool = snakemake.params.get("diags", False)
sv: bool = snakemake.params.get("sv", False)
refine: bool = snakemake.params.get("refine", False)

# Build command line arguments
cmd_args = ["-align", input_seqs, "-output", output_aln]

# Add optional parameters if specified
if maxiters is not None:
    cmd_args.extend(["-maxiters", str(maxiters)])
if diags:
    cmd_args.append("-diags")
if sv:
    cmd_args.append("-sv")
if refine:
    cmd_args.append("-refine")

# Add threads and extra parameters
cmd_args.extend(["-threads", str(threads)])
if extra:
    cmd_args.append(extra)

# Setup logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Execute MUSCLE
shell(f"muscle {' '.join(cmd_args)} {log}") 