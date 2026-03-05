"""Wrapper for CRISPRessoCompare pairwise CRISPR editing comparison tool.

Compares two CRISPResso output directories, producing plots and statistics
summarising differences in editing rates between two conditions.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

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

# --- Inputs: two CRISPResso output directories ---
folder_1 = snakemake.input.folder_1
folder_2 = snakemake.input.folder_2

# --- Resolve output directory and run name ---
out_path = Path(snakemake.output[0])
parent_dir = str(out_path.parent)
bname = out_path.name

name = snakemake.params.get("name", None)
if name is None:
    if bname.startswith("CRISPRessoCompare_"):
        name = bname[len("CRISPRessoCompare_"):]
    else:
        raise ValueError(
            f"Output directory '{bname}' must start with 'CRISPRessoCompare_' "
            "or params.name must be set explicitly."
        )

# --- Optional parameters ---
sample_1_name = snakemake.params.get("sample_1_name", None)
sample_2_name = snakemake.params.get("sample_2_name", None)
reported_qvalue_cutoff = snakemake.params.get("reported_qvalue_cutoff", None)
min_frequency_alleles_around_cut_to_plot = snakemake.params.get(
    "min_frequency_alleles_around_cut_to_plot", None
)
max_rows_alleles_around_cut_to_plot = snakemake.params.get(
    "max_rows_alleles_around_cut_to_plot", None
)
suppress_report = snakemake.params.get("suppress_report", False)
extra = snakemake.params.get("extra", "")

# --- Build command ---
args = [
    f"--crispresso_output_folder_1 {folder_1}",
    f"--crispresso_output_folder_2 {folder_2}",
    f"--output_folder {parent_dir}",
    f"--name {name}",
]

if sample_1_name:
    args.append(f"--sample_1_name {sample_1_name}")
if sample_2_name:
    args.append(f"--sample_2_name {sample_2_name}")
if reported_qvalue_cutoff is not None:
    args.append(f"--reported_qvalue_cutoff {reported_qvalue_cutoff}")
if min_frequency_alleles_around_cut_to_plot is not None:
    args.append(
        f"--min_frequency_alleles_around_cut_to_plot "
        f"{min_frequency_alleles_around_cut_to_plot}"
    )
if max_rows_alleles_around_cut_to_plot is not None:
    args.append(
        f"--max_rows_alleles_around_cut_to_plot "
        f"{max_rows_alleles_around_cut_to_plot}"
    )
if suppress_report:
    args.append("--suppress_report")
if extra:
    args.append(extra)

args_str = " ".join(args)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(f"CRISPRessoCompare {args_str} {log}")
