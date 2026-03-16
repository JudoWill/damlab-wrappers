"""Wrapper for splitting BAM/SAM files by various criteria"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import csv
import os
import re
import sys
import yaml
from collections import Counter
from typing import Dict, Iterator, List, Optional, Tuple

import cigarmath as cm  # type: ignore
import pysam  # type: ignore

if "snakemake" not in locals():
    import snakemake  # type: ignore


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_region(region_str: str) -> Tuple[str, int, int]:
    """Parse 'chr:start-stop' into (chrom, start, end). Coordinates are 0-based half-open."""
    match = re.match(r"^([^:]+):(\d+)-(\d+)$", region_str)
    if not match:
        raise ValueError(
            f"Invalid region format: {region_str!r}. Expected 'chr:start-stop'."
        )
    return match.group(1), int(match.group(2)), int(match.group(3))


def get_sort_order(bam_path: str) -> str:
    """Return the SO field from the BAM header, or 'unknown'."""
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        return bam.header.to_dict().get("HD", {}).get("SO", "unknown")


def build_output_header(bam_path: str) -> pysam.AlignmentHeader:
    """Copy the input BAM header and append a PG line for this wrapper."""
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        header_dict = bam.header.to_dict()
    pg: Dict = {
        "ID": "cigarmath-split",
        "PN": "cigarmath-split",
        "VN": __version__,
        "CL": " ".join(sys.argv),
    }
    if header_dict.get("PG"):
        pg["PP"] = header_dict["PG"][-1]["ID"]
    header_dict.setdefault("PG", []).append(pg)
    return pysam.AlignmentHeader.from_dict(header_dict)


def load_bed_regions(bed_path: str) -> List[Tuple[str, int, int, str]]:
    """Parse a BED file into (chrom, start, end, name) tuples.
    Name is taken from column 4; falls back to 'region_N' if absent.
    """
    regions = []
    with open(bed_path) as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            cols = line.split("\t")
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            name = cols[3].strip() if len(cols) >= 4 else f"region_{i}"
            regions.append((chrom, start, end, name))
    return regions


def load_query_csv(csv_path: str) -> Dict[str, str]:
    """Load a CSV whose first column is query_name and second is category."""
    mapping: Dict[str, str] = {}
    with open(csv_path) as f:
        reader = csv.reader(f)
        header = next(reader, None)  # skip header row
        if header is None:
            return mapping
        for row in reader:
            if len(row) >= 2:
                mapping[row[0].strip()] = row[1].strip()
    return mapping


# ---------------------------------------------------------------------------
# Segment stream helpers
# ---------------------------------------------------------------------------

def iter_units(
    bam_path: str,
    is_name_sorted: bool,
) -> Iterator[Tuple[int, Optional[list], List]]:
    """Yield (ref_start, cigartuples, [segments]) units.

    For name-sorted BAMs the combined_segment_stream merges supplementary
    alignments into a single logical unit.  For coordinate-sorted (or
    unknown-sorted) BAMs each segment is yielded individually.
    """
    stream = cm.io.segment_stream_pysam(bam_path, mode="rb")
    if is_name_sorted:
        valid = (s for s in stream if s.cigartuples)
        for start, cigars, segments in cm.io.combined_segment_stream(valid):
            yield start, cigars, segments
    else:
        for seg in stream:
            yield seg.reference_start, seg.cigartuples, [seg]


# ---------------------------------------------------------------------------
# Required reference overlap filter
# ---------------------------------------------------------------------------

def fully_covers_region(
    start: int,
    cigartuples,
    region_chrom: str,
    region_start: int,
    region_end: int,
    seg_chrom: Optional[str],
) -> bool:
    """Return True if any single mapping block fully spans (region_start, region_end)."""
    if seg_chrom != region_chrom:
        return False
    if cigartuples is None:
        return False
    region_size = region_end - region_start
    for block in cm.reference_mapping_blocks(cigartuples, reference_start=start):
        if cm.block_overlap_length(block, (region_start, region_end)) == region_size:
            return True
    return False


def passes_required_overlap(
    start: int,
    cigartuples,
    segments: List,
    required_regions: List[Tuple[str, int, int]],
) -> bool:
    """Return True only if the read covers ALL required regions."""
    if not required_regions:
        return True
    seg_chrom = segments[0].reference_name
    return all(
        fully_covers_region(start, cigartuples, chrom, rs, re, seg_chrom)
        for chrom, rs, re in required_regions
    )


# ---------------------------------------------------------------------------
# Split criterion functions
# ---------------------------------------------------------------------------

def get_category_flag(
    segments: List,
    flag_categories: Dict[str, Dict[str, bool]],
) -> Optional[str]:
    """Return the first category whose flag conditions all match segments[0]."""
    seg = segments[0]
    for cat_name, conditions in flag_categories.items():
        if all(getattr(seg, attr, None) == val for attr, val in conditions.items()):
            return cat_name
    return None


def get_category_tag(
    segments: List,
    tag_name: str,
    tag_categories: Optional[Dict],
) -> Optional[str]:
    """Return the category for the tag value on segments[0]."""
    try:
        raw = segments[0].get_tag(tag_name)
    except KeyError:
        return None
    if tag_categories:
        return tag_categories.get(str(raw)) or tag_categories.get(raw)
    return str(raw)


def get_category_bed(
    start: int,
    cigartuples,
    segments: List,
    bed_regions: List[Tuple[str, int, int, str]],
) -> Optional[str]:
    """Return the name of the first BED region that the read overlaps."""
    if cigartuples is None:
        return None
    seg_chrom = segments[0].reference_name
    for chrom, reg_start, reg_end, name in bed_regions:
        if seg_chrom != chrom:
            continue
        for block in cm.reference_mapping_blocks(cigartuples, reference_start=start):
            if cm.block_overlap_length(block, (reg_start, reg_end)) > 0:
                return name
    return None


def get_category_size(
    start: int,
    cigartuples,
    size_bins: Dict[str, Tuple],
) -> Optional[str]:
    """Return the first bin name where lo <= query_block_size < hi."""
    if cigartuples is None:
        return None
    q_start, q_end = cm.query_block(cigartuples)
    size = q_end - q_start
    for cat_name, (lo, hi) in size_bins.items():
        lo_val = lo if lo is not None else 0
        if size >= lo_val and (hi is None or size < hi):
            return cat_name
    return None


def get_category_csv(
    query_name: str,
    query_name_map: Dict[str, str],
) -> Optional[str]:
    """Return the category for the query name from a pre-loaded CSV mapping."""
    return query_name_map.get(query_name)


# ---------------------------------------------------------------------------
# Output file management
# ---------------------------------------------------------------------------

def make_output_handle(category: str, output_dir: str, output_format: str, header):
    """Open (or create) an output file handle for the given category."""
    ext = {"fasta": "fa", "fastq": "fq", "sam": "sam", "bam": "bam"}[output_format]
    path = os.path.join(output_dir, f"{category}.{ext}")
    if output_format == "bam":
        return pysam.AlignmentFile(path, "wb", header=header)
    if output_format == "sam":
        return pysam.AlignmentFile(path, "wh", header=header)
    return open(path, "w")


def write_record(handle, segments: List, output_format: str) -> None:
    """Write a read group to the appropriate output handle."""
    if output_format in ("bam", "sam"):
        for seg in segments:
            handle.write(seg)
    else:
        seg = segments[0]
        if not seg.query_sequence:
            return
        name = seg.query_name
        seq = seg.query_sequence
        if output_format == "fasta":
            handle.write(f">{name}\n{seq}\n")
        else:
            quals = seg.query_qualities
            quality = (
                "".join(chr(q + 33) for q in quals)
                if quals is not None
                else "!" * len(seq)
            )
            handle.write(f"@{name}\n{seq}\n+\n{quality}\n")


def _wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i: i + width] for i in range(0, len(seq), width))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# --- Inputs / outputs -------------------------------------------------------
input_bam: str = snakemake.input[0]
output_dir: str = snakemake.output[0]
output_summary_yaml: Optional[str] = (
    snakemake.output[1] if len(snakemake.output) > 1 else None
)

# Optional file inputs
bed_file: Optional[str] = snakemake.input.get("bed_file", None)
csv_file: Optional[str] = snakemake.input.get("csv_file", None)

# --- Parameters -------------------------------------------------------------
split_by: str = snakemake.params.get("split_by")
if not split_by:
    raise ValueError("params.split_by is required.")

output_format: str = snakemake.params.get("output_format", "").lower()
if output_format not in ("fasta", "fastq", "sam", "bam"):
    raise ValueError(
        f"params.output_format must be one of fasta/fastq/sam/bam, got {output_format!r}."
    )

flag_categories: Dict = snakemake.params.get("flag_categories", {})
tag_name: Optional[str] = snakemake.params.get("tag_name", None)
tag_categories: Optional[Dict] = snakemake.params.get("tag_categories", None)
size_bins: Dict = snakemake.params.get("size_bins", {})
top_n: Optional[int] = snakemake.params.get("top_n", None)
sample_name: str = snakemake.params.get("sample_name", "sample")
unclassified_name: str = snakemake.params.get("unclassified_name", "unclassified")
required_overlap_strs: List[str] = snakemake.params.get("required_reference_overlap", []) or []

# Parse required overlap regions
required_regions: List[Tuple[str, int, int]] = [
    parse_region(r) for r in required_overlap_strs
]

# Load auxiliary inputs
bed_regions: List[Tuple[str, int, int, str]] = []
if split_by == "bed":
    if not bed_file:
        raise ValueError("split_by='bed' requires input.bed_file.")
    bed_regions = load_bed_regions(bed_file)

query_name_map: Dict[str, str] = {}
if split_by == "csv":
    if not csv_file:
        raise ValueError("split_by='csv' requires input.csv_file.")
    query_name_map = load_query_csv(csv_file)

# --- Sort order detection ---------------------------------------------------
is_name_sorted: bool = get_sort_order(input_bam) in ("queryname", "collated")

# --- Output header (SAM/BAM only) ------------------------------------------
output_header = None
if output_format in ("sam", "bam"):
    output_header = build_output_header(input_bam)

# --- Create output directory ------------------------------------------------
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Category resolution helper
# ---------------------------------------------------------------------------

def resolve_category(start, cigartuples, segments) -> Optional[str]:
    """Return the split category for this read unit, or None if unresolvable."""
    query_name = segments[0].query_name
    if split_by == "flag":
        return get_category_flag(segments, flag_categories)
    elif split_by == "tag":
        return get_category_tag(segments, tag_name, tag_categories)
    elif split_by == "bed":
        return get_category_bed(start, cigartuples, segments, bed_regions)
    elif split_by == "query_size":
        return get_category_size(start, cigartuples, size_bins)
    elif split_by == "csv":
        return get_category_csv(query_name, query_name_map)
    else:
        raise ValueError(f"Unknown split_by: {split_by!r}")


# ---------------------------------------------------------------------------
# Pass 1 — count categories (only when top_n is requested)
# ---------------------------------------------------------------------------

allowed_categories: Optional[set] = None

if top_n is not None:
    category_counts: Counter = Counter()
    for start, cigartuples, segments in iter_units(input_bam, is_name_sorted):
        if not passes_required_overlap(start, cigartuples, segments, required_regions):
            continue
        cat = resolve_category(start, cigartuples, segments)
        if cat is not None:
            category_counts[cat] += 1
    allowed_categories = {cat for cat, _ in category_counts.most_common(top_n)}

# ---------------------------------------------------------------------------
# Pass 2 (or single pass) — split reads into output files
# ---------------------------------------------------------------------------

output_handles: Dict = {}

total_reads = 0
required_overlap_filtered = 0
category_tally: Counter = Counter()
unclassified_count = 0

for start, cigartuples, segments in iter_units(input_bam, is_name_sorted):
    total_reads += 1

    if not passes_required_overlap(start, cigartuples, segments, required_regions):
        required_overlap_filtered += 1
        continue

    cat = resolve_category(start, cigartuples, segments)

    if cat is None or (allowed_categories is not None and cat not in allowed_categories):
        cat = unclassified_name
        unclassified_count += 1
    else:
        category_tally[cat] += 1

    if cat not in output_handles:
        output_handles[cat] = make_output_handle(cat, output_dir, output_format, output_header)
    write_record(output_handles[cat], segments, output_format)

# Close all output handles
for handle in output_handles.values():
    handle.close()

# ---------------------------------------------------------------------------
# Summary YAML
# ---------------------------------------------------------------------------

if output_summary_yaml:
    summary = {
        "sample_name": sample_name,
        "split_by": split_by,
        "output_format": output_format,
        "output_directory": output_dir,
        "is_name_sorted": is_name_sorted,
        "top_n": top_n,
        "required_reference_overlap": required_overlap_strs or None,
        "total_reads": total_reads,
        "required_overlap_filtered": required_overlap_filtered,
        "categories": dict(category_tally),
        "unclassified": unclassified_count,
    }
    os.makedirs(os.path.dirname(output_summary_yaml) or ".", exist_ok=True)
    with open(output_summary_yaml, "w") as f:
        f.write("# Cigarmath Split\n")
        yaml.dump(summary, f, default_flow_style=False)
