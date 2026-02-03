"""Wrapper for phytreeviz tree visualization.

This wrapper provides a Snakemake interface to phytreeviz, a tool for creating
publication-quality visualizations of phylogenetic trees with customizable
styling options.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.1"

import os
from typing import Optional, Union

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
input_tree: str = snakemake.input[0]
output_plot: str = snakemake.output[0]

# Validate input file exists
if not os.path.exists(input_tree):
    raise FileNotFoundError(f"Input file {input_tree} does not exist")

# Get parameters with defaults
# This is a common pattern in wrappers - providing sensible defaults
extra: str = snakemake.params.get("extra", "")
format: str = snakemake.params.get("format", "newick")  # Default to newick format
width: Union[int, float] = snakemake.params.get("width", 10)  # Default width in inches
height: Union[int, float] = snakemake.params.get("height", 10)  # Default height in inches
show_branch_support: bool = snakemake.params.get("show_branch_support", False)
color_scheme: Optional[str] = snakemake.params.get("color_scheme", None)
dpi: int = snakemake.params.get("dpi", 300)  # Default DPI for raster formats

# Build command line arguments
cmd_args = ["-i", input_tree, "-o", output_plot]

# Add visualization parameters
cmd_args.extend(["--format", format])
cmd_args.extend(["--fig_width", str(width)])
cmd_args.extend(["--fig_height", str(height)])
cmd_args.extend(["--dpi", str(dpi)])

if show_branch_support:
    cmd_args.append("--show_branch_support")
if color_scheme:
    cmd_args.extend(["--color_scheme", color_scheme])

# Add extra parameters if specified
if extra:
    cmd_args.append(extra)

# Setup logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Execute phytreeviz
shell(f"phytreeviz {' '.join(cmd_args)} {log}") 