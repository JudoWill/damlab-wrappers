"""Wrapper for calculating deletion frequencies from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.1.0"

import yaml
import csv
import cigarmath as cm # type: ignore

if "snakemake" not in locals():
    import snakemake # type: ignore

# Get input/output files
input_bam = snakemake.input[0]
output_yaml = snakemake.output[0]
# Check if CSV output is specified
output_csv = snakemake.output[1] if len(snakemake.output) > 1 else None

# Get required parameters
required_region = snakemake.params.get("required_region")
deletion_region = snakemake.params.get("deletion_region")
region_name = snakemake.params.get("region_name", "region")
sample_name = snakemake.params.get("sample_name", "sample")
min_deletion_size = snakemake.params.get("min_deletion_size", 1)

# Validate parameters
if not required_region or not deletion_region:
    raise ValueError("Both required_region and deletion_region parameters must be specified")

# Parse region strings (format: "chr:start-end")
try:
    req_chrom, req_coords = required_region.split(":")
    req_start, req_end = map(int, req_coords.split("-"))
    
    del_chrom, del_coords = deletion_region.split(":")
    del_start, del_end = map(int, del_coords.split("-"))
except ValueError:
    raise ValueError("Region format must be 'chr:start-end'")

segment_stream = cm.io.segment_stream_pysam(input_bam, mode='rb')
fixed_segment_stream = (segment for segment in segment_stream if segment.cigartuples)
combined_segment_stream = cm.io.combined_segment_stream(fixed_segment_stream)

stats = {
    "sample_name": sample_name,
    "region_name": region_name,
    "total_reads": 0,
    "reads_covering_required": 0,
    "reads_covered_partial_deleted": 0,
    "reads_covered_full_deleted": 0,
    "partial_deletion_frequency": None,
    "full_deletion_frequency": None
}

# List to store read-level information if CSV output is requested
read_data = []

req_block_size = req_end - req_start
del_block_size = del_end - del_start
for start, cigars, segments in combined_segment_stream:
    if not cigars:
        continue
    
    stats["total_reads"] += 1
    
    partial_deleted = False
    full_deleted = False

    # Get read name from the first segment
    read_name = segments[0].query_name if segments else "unknown"
    
    # Check if read covers required region
    ref_block = cm.reference_block(cigars, start)
    covers_required = cm.block_overlap_length(ref_block, (req_start, req_end)) == req_block_size
    
    # Only check for deletions if the read covers the required region
    if covers_required:
        stats["reads_covering_required"] += 1
        # Check if read has deletion in specified region
        
        for del_block in cm.reference_deletion_blocks(cigars, start, min_size=min_deletion_size):
            if del_block[1] < del_start:
                # Deletion is before the start of the region, so we can skip
                continue
            if del_block[0] > del_end:
                # All deletions are past the end of the region, so we can break
                break
            if (del_block[0] <= del_start) and (del_block[1] >= del_end):
                stats["reads_covered_full_deleted"] += 1
                full_deleted = True
                break
            if cm.block_overlap_length(del_block, (del_start, del_end)) > min_deletion_size:
                stats["reads_covered_partial_deleted"] += 1
                partial_deleted = True
                break
    
    
    # Store read-level information if CSV output is requested
    if output_csv:
        read_data.append({
            "read_name": read_name,
            "covers_required": covers_required,
            "full_deleted": full_deleted,
            "partial_deleted": partial_deleted
        })

try:
    stats["partial_deletion_frequency"] = stats["reads_covered_partial_deleted"] / stats["reads_covering_required"]
    stats["full_deletion_frequency"] = stats["reads_covered_full_deleted"] / stats["reads_covering_required"]
except ZeroDivisionError:
    stats["partial_deletion_frequency"] = None
    stats["full_deletion_frequency"] = None

# Write YAML output
with open(output_yaml, 'w') as f:
    f.write('# Cigarmath Deletion Frequency\n')
    yaml.dump(stats, f)

# Write CSV output if requested
if output_csv:
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["read_name", "covers_required", "full_deleted", "partial_deleted"])
        writer.writeheader()
        writer.writerows(read_data) 