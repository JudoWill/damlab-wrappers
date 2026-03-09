"""Wrapper for CRISPRessoAggregate multi-run aggregation tool.

Aggregates output from any number of CRISPResso runs into a single combined
report. Each input directory is passed as a separate --prefix argument.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import os
from pathlib import Path

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

# --- Resolve output directory and run name ---
# CRISPRessoAggregate writes CRISPRessoAggregate_on_{name}/ into the cwd,
# so we cd to the parent and derive --name from the output basename.
# Absolute paths are required throughout: after the cd, any relative paths
# (inputs, log redirects) would be resolved against the wrong directory.
out_path = Path(snakemake.output[0]).resolve()
parent_dir = str(out_path.parent)
bname = out_path.name

name = snakemake.params.get("name", None)
if name is None:
    if bname.startswith("CRISPRessoAggregate_on_"):
        name = bname[len("CRISPRessoAggregate_on_"):]
    else:
        raise ValueError(
            f"Output directory '{bname}' must start with 'CRISPRessoAggregate_on_' "
            "or params.name must be set explicitly."
        )

# --- Inputs: one or more CRISPResso output directories ---
# Accept either a named list (input.crispresso_dirs) or positional inputs.
if hasattr(snakemake.input, "crispresso_dirs"):
    input_dirs = snakemake.input.crispresso_dirs
    if isinstance(input_dirs, str):
        input_dirs = [input_dirs]
else:
    input_dirs = list(snakemake.input)

if not input_dirs:
    raise ValueError("At least one CRISPResso output directory must be provided as input.")

# Resolve to absolute paths so the cd to parent_dir doesn't invalidate them
input_dirs = [os.path.abspath(d) for d in input_dirs]

# --- Optional parameters ---
min_reads_for_inclusion = snakemake.params.get("min_reads_for_inclusion", None)
max_samples_per_summary_plot = snakemake.params.get("max_samples_per_summary_plot", None)
suppress_report = snakemake.params.get("suppress_report", False)
suppress_plots = snakemake.params.get("suppress_plots", False)
threads = snakemake.threads
extra = snakemake.params.get("extra", "")

# --- Build command ---
# Each input directory is passed as a separate --prefix so CRISPRessoAggregate
# globs exactly that directory path.
prefix_args = " ".join(f"--prefix {d}" for d in input_dirs)

args = [
    prefix_args,
    f"--name {name}",
    f"--n_processes {threads}",
]

if min_reads_for_inclusion is not None:
    args.append(f"--min_reads_for_inclusion {min_reads_for_inclusion}")
if max_samples_per_summary_plot is not None:
    args.append(f"--max_samples_per_summary_plot {max_samples_per_summary_plot}")
if suppress_report:
    args.append("--suppress_report")
if suppress_plots:
    args.append("--suppress_plots")
if extra:
    args.append(extra)

args_str = " ".join(args)

# log_fmt_shell may return relative paths; make them absolute so the cd
# to parent_dir doesn't cause bash to look for the log in the wrong place.
if snakemake.log:
    log_path = os.path.abspath(str(snakemake.log[0]))
    log = f"> {log_path} 2>&1"
else:
    log = ""

# cd to parent so CRISPRessoAggregate writes its output folder there
shell(f"cd {parent_dir} && CRISPRessoAggregate {args_str} {log}")
