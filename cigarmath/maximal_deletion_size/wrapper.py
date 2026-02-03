"""Wrapper for calculating maximal deletion sizes from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import yaml
import csv
import cigarmath as cm # type: ignore

if "snakemake" not in locals():
    import snakemake # type: ignore

# Get input/output files
input_bam = snakemake.input[0]
output_yaml = snakemake.output[0]
output_csv = snakemake.output[1]

# Get parameters
sample_name = snakemake.params.get("sample_name", "sample")
chrom = snakemake.params.get("chrom")
min_size = snakemake.params.get("min_size", 1)
all_query_names = set()

def _capture_query_names(stream):
    """Capture all query names from the stream"""

    for item in stream:
        all_query_names.add(item.query_name)
        yield item


segment_stream = _capture_query_names(cm.io.segment_stream_pysam(input_bam, mode='rb'))
combined_segment_stream = cm.io.combined_segment_stream(segment_stream)





# Initialize metrics
metrics = {
    "sample_name": sample_name,
    "mapped_reads": 0,
    "average_ref_block_size": 0.0,
    "average_largest_deletion_size": 0.0
}

# List to store read-level information
read_data = []

total_ref_block_size = 0
total_largest_deletion_size = 0

for start, cigars, segments in combined_segment_stream:
    if not cigars:
        continue
    
    # Get read name from the first segment
    read_name = segments[0].query_name if segments else "unknown"
    
    # Skip if chromosome filter is specified and read doesn't match
    if chrom and segments[0].reference_name != chrom:
        continue
    
    # Calculate reference block
    ref_block = cm.reference_block(cigars, start)
    ref_block_size = ref_block[1] - ref_block[0]
    
    # Find largest deletion
    largest_deletion = None
    largest_deletion_size = 0
    
    for del_block in cm.reference_deletion_blocks(cigars, reference_start=start, min_size=min_size):
        del_size = del_block[1] - del_block[0]
        if del_size > largest_deletion_size:
            largest_deletion = del_block
            largest_deletion_size = del_size
    
    # Update metrics
    if ref_block_size > 0:  # Only count mapped reads
        metrics["mapped_reads"] += 1
        total_ref_block_size += ref_block_size
        total_largest_deletion_size += largest_deletion_size
    
    # Store read-level information
    read_data.append({
        "read_name": read_name,
        "reference_block_start": ref_block[0],
        "reference_block_stop": ref_block[1],
        "reference_block_size": ref_block_size,
        "largest_deletion_start": largest_deletion[0] if largest_deletion else None,
        "largest_deletion_stop": largest_deletion[1] if largest_deletion else None,
        "largest_deletion_size": largest_deletion_size
    })

metrics["total_reads"] = len(all_query_names)

# Calculate averages
if metrics["mapped_reads"] > 0:
    metrics["average_ref_block_size"] = total_ref_block_size / metrics["mapped_reads"]
    metrics["average_largest_deletion_size"] = total_largest_deletion_size / metrics["mapped_reads"]

# Write YAML output
with open(output_yaml, 'w') as f:
    f.write('# Cigarmath Maximal Deletion Size\n')
    yaml.dump(metrics, f)

# Write CSV output
with open(output_csv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=[
        "read_name", "reference_block_start", "reference_block_stop", 
        "reference_block_size", "largest_deletion_start", "largest_deletion_stop", 
        "largest_deletion_size"
    ])
    writer.writeheader()
    writer.writerows(read_data) 