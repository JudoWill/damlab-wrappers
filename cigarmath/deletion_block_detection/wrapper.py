"""Wrapper for detecting deletion blocks from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.1.0"

import math
import yaml
import csv
from collections import Counter, defaultdict
from typing import Dict, Set, Tuple, List, Optional
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


def merge_blocks(
    blocks: Counter,
    n: int,
    *,
    sum_counts: bool = True,
) -> Tuple[Counter, Dict[Tuple[int, int], Tuple[int, int]]]:
    """
    Merge blocks whose start and stop are within n units of each other.
    
    Args:
        blocks: Counter keyed by (start, stop) tuples.
        n: integer proximity threshold.
        sum_counts: if True, representative's count is sum of component counts.
                    if False, representative's count is the original count of the chosen block.
    
    Returns:
        Tuple of:
          - Counter of merged blocks: representative_block -> count
          - Dict mapping each original block to its representative block
    """
    if not blocks:
        return Counter(), {}

    items: List[Tuple[Tuple[int, int], int]] = list(blocks.items())
    m = len(items)

    # Union-find
    parent = list(range(m))
    rank = [0] * m

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[rb] < rank[ra]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    for i in range(m):
        (si, ei), _ = items[i]
        for j in range(i + 1, m):
            (sj, ej), _ = items[j]
            if abs(si - sj) <= n and abs(ei - ej) <= n:
                union(i, j)

    groups: Dict[int, List[int]] = defaultdict(list)
    for idx in range(m):
        groups[find(idx)].append(idx)

    merged: Counter = Counter()
    remap: Dict[Tuple[int, int], Tuple[int, int]] = {}
    for group_idxs in groups.values():
        rep_idx = max(
            group_idxs,
            key=lambda idx: (items[idx][1], -items[idx][0][0], -items[idx][0][1]),
        )
        rep_key = items[rep_idx][0]
        if sum_counts:
            total = sum(items[idx][1] for idx in group_idxs)
            merged[rep_key] += total
        else:
            merged[rep_key] += items[rep_idx][1]
        for idx in group_idxs:
            remap[items[idx][0]] = rep_key

    return merged, remap


def calculate_shannon_entropy(counter: Counter) -> float:
    """Calculate Shannon entropy of a Counter's value distribution."""
    total = sum(counter.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counter.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log(p)
    return entropy


# Get input/output files
input_bams = snakemake.input.bams if hasattr(snakemake.input, 'bams') else list(snakemake.input)
output_reads_csv = snakemake.output.reads
output_deletions_csv = snakemake.output.deletions
output_summary_yaml = snakemake.output.summary

# Get parameters
min_deletion_size = snakemake.params.get("min_deletion_size", 50)
merge_distance = snakemake.params.get("merge_distance", 0)
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

# Optionally merge deletion blocks with similar endpoints
if merge_distance > 0:
    deletion_counter, block_remap = merge_blocks(deletion_counter, merge_distance)

    # Update deletion coordinates in read_data to use representative blocks
    for read in read_data:
        if read['deletions']:
            old_dels: List[Tuple[int, int]] = [
                tuple(map(int, d.split('-'))) for d in read['deletions'].split(';')  # type: ignore
            ]
            new_dels = [block_remap.get(d, d) for d in old_dels]
            # Deduplicate while preserving order (two old blocks may collapse to the same rep)
            seen_dels: Set[Tuple[int, int]] = set()
            unique_dels = []
            for d in new_dels:
                if d not in seen_dels:
                    seen_dels.add(d)
                    unique_dels.append(d)
            read['deletions'] = ';'.join(f"{d[0]}-{d[1]}" for d in unique_dels)

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

# Calculate summary statistics
total_deletion_count = sum(deletion_counter.values())
richness = len(deletion_counter)
shannon_entropy = calculate_shannon_entropy(deletion_counter)

# Write summary YAML for MultiQC
summary = {
    'sample_name': sample_name,
    'total_reads': total_reads,
    'reads_with_deletions': reads_with_deletions,
    'unique_deletion_count': richness,
    'total_deletion_count': total_deletion_count,
    'deletion_frequency': reads_with_deletions / total_reads if total_reads > 0 else 0.0,
    'deletion_richness': richness,
    'deletion_shannon_entropy': round(shannon_entropy, 6),
    'min_deletion_size': min_deletion_size,
    'merge_distance': merge_distance,
    'input_bam_count': len(input_bams),
    'allowedlist_used': allowed_deletions is not None,
    'allowedlist_size': len(allowed_deletions) if allowed_deletions else 0
}

with open(output_summary_yaml, 'w') as f:
    f.write('# Cigarmath Deletion Block Detection\n')
    yaml.dump(summary, f, default_flow_style=False)
