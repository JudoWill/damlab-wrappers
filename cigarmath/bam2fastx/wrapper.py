"""Wrapper for converting BAM/SAM files to FASTA or FASTQ"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "0.0.1"

import os
from typing import Optional, Iterator

import cigarmath as cm

if "snakemake" not in locals():
    import snakemake


def check_bam_index(path: str) -> bool:
    """Check if BAM index (.bai) file exists."""
    if not path.lower().endswith(".bam"):
        return False
    index_path = f"{path}.bai"
    alt_index_path = f"{os.path.splitext(path)[0]}.bai"
    return os.path.exists(index_path) or os.path.exists(alt_index_path)


def get_segment_stream(
    path: str,
    region: Optional[str],
    min_mapq: int = 0,
) -> Iterator:
    """Get a stream of aligned segments, using fetch if BAM index exists."""
    has_index = check_bam_index(path)
    use_fetch = region if (has_index and region) else None

    return cm.io.segment_stream_pysam(
        path,
        mode="rb" if path.lower().endswith(".bam") else "r",
        fetch=use_fetch,
        min_mapq=min_mapq,
    )


def _infer_output_format(output_path: str, param_format: Optional[str]) -> str:
    """Infer output format from file extension or explicit param."""
    if param_format:
        fmt = param_format.lower()
        if fmt in ("fasta", "fa"):
            return "fasta"
        if fmt in ("fastq", "fq"):
            return "fastq"
        raise ValueError(f"output_format must be 'fasta' or 'fastq', got {param_format!r}")

    ext = os.path.splitext(output_path)[1].lower()
    if ext in (".fa", ".fasta"):
        return "fasta"
    if ext in (".fq", ".fastq"):
        return "fastq"
    raise ValueError(
        f"Cannot infer output format from extension {ext!r}. "
        "Use .fa/.fasta for FASTA or .fq/.fastq for FASTQ, or set params.output_format."
    )


def _wrap_fasta(seq: str, width: int = 60) -> str:
    """Wrap sequence to width characters per line."""
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


# Get input/output files
input_path = snakemake.input[0]
output_path = snakemake.output[0]

# Get parameters
mapped_only = snakemake.params.get("mapped_only", False)
region = snakemake.params.get("region")
min_mapq = snakemake.params.get("min_mapq", 0)
output_format_param = snakemake.params.get("output_format")

output_format = _infer_output_format(output_path, output_format_param)

# Stream segments
segment_stream = get_segment_stream(input_path, region, min_mapq)

with open(output_path, "w") as f:
    for segment in segment_stream:
        if mapped_only and segment.is_unmapped:
            continue
        if not segment.query_sequence:
            continue

        name = segment.query_name
        seq = segment.query_sequence

        if output_format == "fasta":
            f.write(f">{name}\n")
            f.write(_wrap_fasta(seq))
            f.write("\n")
        else:
            # FASTQ
            quals = segment.query_qualities
            if quals is not None:
                quality = "".join(chr(q + 33) for q in quals)
            else:
                quality = "!" * len(seq)
            f.write(f"@{name}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{quality}\n")
