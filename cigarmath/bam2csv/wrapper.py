"""Wrapper for converting BAM/SAM files to CSV with specified fields"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import pysam
import csv
from typing import List, Dict, Any

if "snakemake" not in locals():
    import snakemake

def get_field_value(segment: pysam.AlignedSegment, field: str) -> Any:
    """Get the value of a field from a BAM/SAM segment.
    
    Args:
        segment: pysam AlignedSegment object
        field: Field name to extract (e.g., 'query_name', 'query_sequence', or a tag)
        
    Returns:
        The value of the requested field
    """
    # Handle standard BAM fields
    if hasattr(segment, field):
        return getattr(segment, field)
    
    # Handle tags (e.g., 'NM:i:1' -> 'NM')
    if field in segment.tags:
        return segment.get_tag(field)
    
    return None

# Get input/output files
input_bam = snakemake.input[0]
output_csv = snakemake.output[0]

# Get required parameters
fields = snakemake.params.get("fields", ["query_name", "reference_name", "reference_start"])
if not isinstance(fields, list):
    raise ValueError("params['fields'] must be a list of field names")

# Open BAM file
with pysam.AlignmentFile(input_bam) as bam:
    # Write output CSV
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow(fields)
        
        # Process each segment
        for segment in bam:
            # Get values for each field
            row = []
            for field in fields:
                value = get_field_value(segment, field)
                row.append(value)
            
            # Write row
            writer.writerow(row) 