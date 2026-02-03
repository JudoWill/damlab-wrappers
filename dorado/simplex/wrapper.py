"""Wrapper for dorado simplex basecalling"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

from snakemake.shell import shell # type: ignore
from os import path

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact function
    import snakemake # type: ignore

# Extract arguments from snakemake object
pod5_file = snakemake.input.pod
output_bam = snakemake.output[0]

# Get optional parameters
model = snakemake.params.get("model", "sup")
models_dir = snakemake.params.get("models_directory", "")
dorado_path = snakemake.params.get("dorado_path", "dorado")  # default to 'dorado' in PATH
models_dir_arg = f"--models-directory {models_dir}" if models_dir else ""

# Handle recursive mode for directories
recursive_arg = ""
if path.isdir(pod5_file):
    recursive_arg = "--recursive"

# Handle reference - now optional
reference_arg = ""
if hasattr(snakemake.input, 'reference'):
    reference_arg = f"--reference {snakemake.input.reference}"

# Handle GPU resources
gpu_arg = "--device cuda:all"  # default to all GPUs
if 'gpu' in snakemake.params.keys():
    gpu = snakemake.params.get('gpu', 'all')
    if gpu == 'all':
        gpu_arg = "--device cuda:all"
    else:
        gpu_arg = f"--device cuda:{gpu}"


# Handle logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Build and execute command
shell(
    f"{dorado_path} basecaller"
    f" {gpu_arg}"
    f" {models_dir_arg}"
    f" {reference_arg}"
    f" {recursive_arg}"
    f" {model}"
    f" {pod5_file}"
    f" > {output_bam}"
    f" {log}"
) 