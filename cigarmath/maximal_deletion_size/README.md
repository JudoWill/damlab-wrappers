# Maximal Deletion Size Wrapper

This Snakemake wrapper calculates the maximal deletion size for each read in a BAM file, along with reference mapping block information.

## Input

- BAM file (sorted by read name)

## Output

1. YAML file containing summary metrics:
   - Total number of reads
   - Number of mapped reads
   - Average reference block size
   - Average largest deletion size

2. CSV file containing per-read information:
   - Read name
   - Reference block start position
   - Reference block stop position
   - Reference block size
   - Largest deletion start position
   - Largest deletion stop position
   - Largest deletion size

## Parameters

- `min_size` (optional): The smallest deletion size considered for the block operation (default: 1). 
- `sample_name` (optional): Name of the sample for the metrics file (default: "sample")


## Example Usage

```python
rule maximal_deletion_size:
    input:
        bam="path/to/input.bam"
    output:
        yaml="path/to/output.yaml",
        csv="path/to/output.csv"
    params:
        sample_name="sample1"
    wrapper:
        "damlab-wrappers/cigarmath/maximal_deletion_size"
```

## Notes

- The BAM file should be sorted by read name
- The wrapper groups aligned segments of the same query name together
- Reference block size is calculated as the total span of the read's alignment
- The largest deletion is determined by comparing all deletions in the read's CIGAR string 