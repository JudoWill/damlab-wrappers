"""Wrapper for detecting deletion blocks from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import yaml
import csv
from collections import Counter
from typing import Set, Tuple, List, Optional
import cigarmath as cm  # type: ignore

if "snakemake" not in locals():
    import snakemake  # type: ignore


def parse_allowedlist(allowedlist_path: str) -> Set[Tuple[int, int]]:
    """Parse allowedlist file containing deletion coordinates.
    
    Expected format: CSV with columns 'start' and 'end' or 'deletion_start' and 'deletion_end'
    """
    allowed = set()
    with open(allowedlist_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'start' in row and 'end' in row:
                allowed.add((int(row['start']), int(row['end'])))
            elif 'deletion_start' in row and 'deletion_end' in row:
                allowed.add((int(row['deletion_start']), int(row['deletion_end'])))
    return allowed


def get_combined_segments(bam_paths: List[str]):
    """Get combined segments from one or more BAM files."""
    all_segments = []
    for bam_path in bam_paths:
        segments = list(cm.io.segment_stream_pysam(bam_path, mode='rb'))
        all_segments.extend(segments)
    
    all_segments.sort(key=lambda x: x.query_name)
    valid_segments = [s for s in all_segments if s.cigartuples]
    return list(cm.io.combined_segment_stream(iter(valid_segments)))


# Get input/output files
input_bams = snakemake.input.bams if hasattr(snakemake.input, 'bams') else list(snakemake.input)
output_reads_csv = snakemake.output.reads
output_deletions_csv = snakemake.output.deletions
output_summary_yaml = snakemake.output.summary

# Get parameters
min_deletion_size = snakemake.params.get("min_deletion_size", 50)
sample_name = snakemake.params.get("sample_name", "sample")

# Get optional allowedlist
allowedlist_input = snakemake.input.get("allowedlist", None)
allowed_deletions: Optional[Set[Tuple[int, int]]] = None
if allowedlist_input:
    allowed_deletions = parse_allowedlist(allowedlist_input)

# Ensure input_bams is a list
if isinstance(input_bams, str):
    input_bams = [input_bams]

# Process all BAM files together
combined_stream = get_combined_segments(input_bams)

# Data collection
read_data: List[dict] = []
deletion_counter: Counter = Counter()
total_reads = 0
reads_with_deletions = 0

# First pass: collect all deletions and read data
for start, cigartuples, segments in combined_stream:
    if not cigartuples:
        continue
    
    total_reads += 1
    read_name = segments[0].query_name
    ref_end = start + cm.reference_offset(cigartuples)
    
    deletions = list(cm.reference_deletion_blocks(
        cigartuples, 
        reference_start=start, 
        min_size=min_deletion_size
    ))
    
    # Filter by allowedlist if specified
    if allowed_deletions is not None:
        deletions = [d for d in deletions if d in allowed_deletions]
    
    if deletions:
        reads_with_deletions += 1
    
    read_data.append({
        'read_name': read_name,
        'reference_start': start,
        'reference_end': ref_end,
        'deletions': ';'.join(f"{d[0]}-{d[1]}" for d in deletions) if deletions else ''
    })
    
    for del_block in deletions:
        deletion_counter[del_block] += 1

# Second pass: calculate coverage for each deletion
coverage_counter: Counter = Counter()
for read in read_data:
    read_start = read['reference_start']
    read_end = read['reference_end']
    for del_block in deletion_counter:
        if read_start <= del_block[0] and read_end >= del_block[1]:
            coverage_counter[del_block] += 1

# Write read-centered CSV
with open(output_reads_csv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['read_name', 'reference_start', 'reference_end', 'deletions'])
    writer.writeheader()
    writer.writerows(read_data)

# Write deletion-centered CSV
deletion_data = []
for del_block, read_count in deletion_counter.most_common():
    deletion_data.append({
        'deletion_start': del_block[0],
        'deletion_end': del_block[1],
        'deletion_size': del_block[1] - del_block[0],
        'read_count': read_count,
        'coverage_count': coverage_counter[del_block]
    })

with open(output_deletions_csv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['deletion_start', 'deletion_end', 'deletion_size', 'read_count', 'coverage_count'])
    writer.writeheader()
    writer.writerows(deletion_data)

# Calculate total deletion count (sum across all reads)
total_deletion_count = sum(deletion_counter.values())

# Write summary YAML for MultiQC
summary = {
    'sample_name': sample_name,
    'total_reads': total_reads,
    'reads_with_deletions': reads_with_deletions,
    'unique_deletion_count': len(deletion_counter),
    'total_deletion_count': total_deletion_count,
    'deletion_frequency': reads_with_deletions / total_reads if total_reads > 0 else 0.0,
    'min_deletion_size': min_deletion_size,
    'input_bam_count': len(input_bams),
    'allowedlist_used': allowed_deletions is not None,
    'allowedlist_size': len(allowed_deletions) if allowed_deletions else 0
}

with open(output_summary_yaml, 'w') as f:
    f.write('# Cigarmath Deletion Block Detection\n')
    yaml.dump(summary, f, default_flow_style=False)
