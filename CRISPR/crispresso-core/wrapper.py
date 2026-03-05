"""Wrapper for CRISPResso CRISPR editing analysis tool.

Runs CRISPResso on a single amplicon. Amplicon sequence may be provided
either as a string via params.amplicon_seq or as a FASTA file via
input.amplicon_fasta.
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


def read_fasta(path):
    """Return list of (id, seq) tuples parsed from a FASTA file."""
    records = []
    name, seq_parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts)))
                name = line[1:].split()[0]
                seq_parts = []
            elif line:
                seq_parts.append(line)
        if name is not None:
            records.append((name, "".join(seq_parts)))
    return records


# --- Resolve amplicon sequence ---
amplicon_fasta = snakemake.input.get("amplicon_fasta", None)
if amplicon_fasta:
    records = read_fasta(amplicon_fasta)
    if not records:
        raise ValueError(f"No sequences found in amplicon FASTA: {amplicon_fasta}")
    amplicon_seq = ",".join(seq for _, seq in records)
    default_amplicon_name = ",".join(rec_id for rec_id, _ in records)
else:
    amplicon_seq = snakemake.params.get("amplicon_seq", None)
    if not amplicon_seq:
        raise ValueError(
            "Either input.amplicon_fasta or params.amplicon_seq must be provided."
        )
    default_amplicon_name = None

# --- Resolve output directory and run name ---
out_path = Path(snakemake.output[0])
parent_dir = str(out_path.parent)
bname = out_path.name

name = snakemake.params.get("name", None)
if name is None:
    if bname.startswith("CRISPResso_on_"):
        name = bname[len("CRISPResso_on_"):]
    else:
        raise ValueError(
            f"Output directory '{bname}' must start with 'CRISPResso_on_' "
            "or params.name must be set explicitly."
        )

# --- Inputs ---
fastq_r1 = snakemake.input.fastq_r1
fastq_r2 = snakemake.input.get("fastq_r2", None)

# --- Optional parameters ---
amplicon_name = snakemake.params.get("amplicon_name", default_amplicon_name)
guide_seq = snakemake.params.get("guide_seq", None)
guide_name = snakemake.params.get("guide_name", None)
expected_hdr_amplicon_seq = snakemake.params.get("expected_hdr_amplicon_seq", None)
quantification_window_size = snakemake.params.get("quantification_window_size", None)
quantification_window_center = snakemake.params.get("quantification_window_center", None)
min_average_read_quality = snakemake.params.get("min_average_read_quality", None)
min_single_bp_quality = snakemake.params.get("min_single_bp_quality", None)
ignore_substitutions = snakemake.params.get("ignore_substitutions", False)
ignore_insertions = snakemake.params.get("ignore_insertions", False)
ignore_deletions = snakemake.params.get("ignore_deletions", False)
discard_indel_reads = snakemake.params.get("discard_indel_reads", False)
base_editor_output = snakemake.params.get("base_editor_output", False)
conversion_nuc_from = snakemake.params.get("conversion_nuc_from", None)
conversion_nuc_to = snakemake.params.get("conversion_nuc_to", None)
extra = snakemake.params.get("extra", "")
threads = snakemake.threads

# --- Build command ---
args = [
    f"--fastq_r1 {fastq_r1}",
    f"--amplicon_seq {amplicon_seq}",
    f"--output_folder {parent_dir}",
    f"--name {name}",
    f"--n_processes {threads}",
]

if fastq_r2:
    args.append(f"--fastq_r2 {fastq_r2}")
if amplicon_name:
    args.append(f"--amplicon_name {amplicon_name}")
if guide_seq:
    args.append(f"--guide_seq {guide_seq}")
if guide_name:
    args.append(f"--guide_name {guide_name}")
if expected_hdr_amplicon_seq:
    args.append(f"--expected_hdr_amplicon_seq {expected_hdr_amplicon_seq}")
if quantification_window_size is not None:
    args.append(f"--quantification_window_size {quantification_window_size}")
if quantification_window_center is not None:
    args.append(f"--quantification_window_center {quantification_window_center}")
if min_average_read_quality is not None:
    args.append(f"--min_average_read_quality {min_average_read_quality}")
if min_single_bp_quality is not None:
    args.append(f"--min_single_bp_quality {min_single_bp_quality}")
if ignore_substitutions:
    args.append("--ignore_substitutions")
if ignore_insertions:
    args.append("--ignore_insertions")
if ignore_deletions:
    args.append("--ignore_deletions")
if discard_indel_reads:
    args.append("--discard_indel_reads")
if base_editor_output:
    args.append("--base_editor_output")
if conversion_nuc_from:
    args.append(f"--conversion_nuc_from {conversion_nuc_from}")
if conversion_nuc_to:
    args.append(f"--conversion_nuc_to {conversion_nuc_to}")
if extra:
    args.append(extra)

args_str = " ".join(args)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(f"CRISPResso {args_str} {log}")
