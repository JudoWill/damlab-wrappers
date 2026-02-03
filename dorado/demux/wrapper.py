"""Wrapper for dorado demux"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

from snakemake.shell import shell # type: ignore
from os import path
import os
import glob
import shutil
import tempfile

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact funtion
    import snakemake # type: ignore

# Extract arguments from snakemake object
reads = snakemake.input.reads
output_dir = snakemake.output[0]

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Get required parameters
kit_name = snakemake.params.get("kit_name", "")
if not kit_name and not snakemake.params.get("no_classify", False):
    raise ValueError("Either kit_name or no_classify must be specified")

# Get optional parameters
dorado_path = snakemake.params.get("dorado_path", "dorado")
sample_sheet = snakemake.params.get("sample_sheet", "")
max_reads = snakemake.params.get("max_reads", 0)
read_ids = snakemake.params.get("read_ids", "")
emit_fastq = snakemake.params.get("emit_fastq", False)
emit_summary = snakemake.params.get("emit_summary", False)
barcode_both_ends = snakemake.params.get("barcode_both_ends", False)
no_trim = snakemake.params.get("no_trim", False)
sort_bam = snakemake.params.get("sort_bam", False)
barcode_arrangement = snakemake.params.get("barcode_arrangement", "")
barcode_sequences = snakemake.params.get("barcode_sequences", "")
barcode_to_output = snakemake.params.get("barcode_to_output", {})

# Create a temporary directory using tempfile
with tempfile.TemporaryDirectory(dir=snakemake.params.get("tempdir", None)) as temp_output:
    # Handle recursive mode for directories
    recursive_arg = "--recursive" if path.isdir(reads) else ""

    # Build parameter strings
    kit_arg = f"--kit-name {kit_name}" if kit_name else "--no-classify"
    sample_sheet_arg = f"--sample-sheet {sample_sheet}" if sample_sheet else ""
    max_reads_arg = f"--max-reads {max_reads}" if max_reads else ""
    read_ids_arg = f"--read-ids {read_ids}" if read_ids else ""
    emit_fastq_arg = "--emit-fastq" if emit_fastq else ""
    emit_summary_arg = "--emit-summary" if emit_summary else ""
    barcode_both_ends_arg = "--barcode-both-ends" if barcode_both_ends else ""
    no_trim_arg = "--no-trim" if no_trim else ""
    sort_bam_arg = "--sort-bam" if sort_bam else ""
    barcode_arrangement_arg = f"--barcode-arrangement {barcode_arrangement}" if barcode_arrangement else ""
    barcode_sequences_arg = f"--barcode-sequences {barcode_sequences}" if barcode_sequences else ""

    # Get number of threads
    threads_arg = f"--threads {snakemake.threads}" if snakemake.threads > 1 else ""

    # Handle logging
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    # Build and execute command
    shell(
        f"{dorado_path} demux"
        f" {recursive_arg}"
        f" --output-dir {temp_output}"
        f" {kit_arg}"
        f" {sample_sheet_arg}"
        f" {threads_arg}"
        f" {max_reads_arg}"
        f" {read_ids_arg}"
        f" {emit_fastq_arg}"
        f" {emit_summary_arg}"
        f" {barcode_both_ends_arg}"
        f" {no_trim_arg}"
        f" {sort_bam_arg}"
        f" {barcode_arrangement_arg}"
        f" {barcode_sequences_arg}"
        f" {reads}"
        f" {log}"
    )

    # Move and rename files according to barcode mapping
    extension = ".fastq" if emit_fastq else ".bam"
    for file in glob.glob(f"{temp_output}/*{extension}"):
        filename = os.path.basename(file)
        print(filename)
        # Extract barcode from filename
        for barcode, output_name in barcode_to_output.items():
            if f"_{barcode}" in filename:
                print(f"Found {barcode} in {filename}")
                target_file = os.path.join(output_dir, output_name) + extension
                if path.exists(target_file):
                    os.remove(target_file)
                shutil.move(file, target_file)
                barcode_to_output.pop(barcode)
                break
    
    # Anything left in barcode_to_output is missing
    # Create empty files for them to keep snakemake happy
    for barcode, output_name in barcode_to_output.items():
        target_file = os.path.join(output_dir, output_name) + extension
        if path.exists(target_file):
            os.remove(target_file)
        shell(f"touch {target_file}")
