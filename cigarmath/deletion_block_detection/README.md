# Deletion Block Detection

A wrapper for detecting deletion blocks from BAM files using the cigarmath library. This wrapper analyzes aligned reads to identify deletions of a specified minimum size and outputs both read-centered and deletion-centered statistics.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

### Basic Usage

```python
rule detect_deletion_blocks:
    input:
        bams="sorted.bam"
    output:
        reads="deletion_reads.csv",
        deletions="deletion_blocks.csv",
        summary="deletion_summary.yaml"
    params:
        min_deletion_size=50,
        sample_name="patient1"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/deletion_block_detection"
```

### Multiple BAM Files

```python
rule detect_deletion_blocks:
    input:
        bams=["sample1.bam", "sample2.bam", "sample3.bam"]
    output:
        reads="deletion_reads.csv",
        deletions="deletion_blocks.csv",
        summary="deletion_summary.yaml"
    params:
        min_deletion_size=50,
        sample_name="combined_samples"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/deletion_block_detection"
```

### With Allowedlist Filter

```python
rule detect_deletion_blocks:
    input:
        bams="sorted.bam",
        allowedlist="allowed_deletions.csv"
    output:
        reads="deletion_reads.csv",
        deletions="deletion_blocks.csv",
        summary="deletion_summary.yaml"
    params:
        min_deletion_size=50,
        sample_name="filtered_sample"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/deletion_block_detection"
```

## Parameters

- `min_deletion_size` (int, optional): Minimum size of deletions to detect. Defaults to 50.
- `sample_name` (str, optional): Sample name to include in metrics. Defaults to "sample".

## Input

- `bams`: Single BAM file or list of BAM files containing aligned reads
- `allowedlist` (optional): CSV file containing allowed deletions with columns `start`/`end` or `deletion_start`/`deletion_end`

## Output

### Read-centered CSV (`reads`)

Contains one row per read with deletion information:

| Column | Description |
|--------|-------------|
| `read_name` | Name of the read |
| `reference_start` | Start position of alignment on reference |
| `reference_end` | End position of alignment on reference |
| `deletions` | Semicolon-separated list of deletions in format `start-end` |

### Deletion-centered CSV (`deletions`)

Contains one row per unique deletion:

| Column | Description |
|--------|-------------|
| `deletion_start` | Start position of deletion on reference |
| `deletion_end` | End position of deletion on reference |
| `deletion_size` | Size of deletion (end - start) |
| `read_count` | Number of reads containing this deletion |
| `coverage_count` | Number of reads fully covering this deletion region |

### Summary YAML (`summary`)

Contains summary statistics for MultiQC integration:

- `sample_name`: Sample identifier
- `total_reads`: Total number of reads processed
- `reads_with_deletions`: Number of reads containing at least one deletion
- `unique_deletion_count`: Number of unique deletion blocks found
- `total_deletion_count`: Total number of deletions across all reads
- `deletion_frequency`: Fraction of reads with deletions
- `min_deletion_size`: Minimum deletion size threshold used
- `input_bam_count`: Number of input BAM files processed
- `allowedlist_used`: Whether an allowedlist filter was applied
- `allowedlist_size`: Number of deletions in allowedlist (if used)

## Allowedlist Format

The optional allowedlist CSV should have either:
- Columns `start` and `end`, or
- Columns `deletion_start` and `deletion_end`

Example:
```csv
start,end
1084,3594
1075,6368
6903,7666
```

When an allowedlist is provided, only deletions matching entries in the list will be reported.

## Error Handling

The wrapper includes error handling for:
- Invalid BAM files
- Missing required outputs
- Empty or malformed allowedlist files

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors

* Will Dampier, PhD

## Software Requirements

* [pysam](https://pysam.readthedocs.io/) (tested with v0.19.0)
* [cigarmath](https://github.com/DamLabResources/cigarmath) (tested with v0.1.0)
