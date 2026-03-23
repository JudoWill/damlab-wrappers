# Deletion Block Detection

A wrapper for detecting deletion blocks from BAM files using the cigarmath library. This wrapper analyzes aligned reads to identify deletions of a specified minimum size and outputs both read-centered and deletion-centered statistics.

## Version

Current version: 1.2.0

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
        summary="deletion_summary.yaml",
        query_stats="deletion_query_stats.csv",
    params:
        min_deletion_size=50,
        sample_name="patient1"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/deletion_block_detection"
```

`query_stats` is always produced: header-only when `query` is omitted; see [Query regions](#query-regions).

### Multiple BAM Files

```python
rule detect_deletion_blocks:
    input:
        bams=["sample1.bam", "sample2.bam", "sample3.bam"]
    output:
        reads="deletion_reads.csv",
        deletions="deletion_blocks.csv",
        summary="deletion_summary.yaml",
        query_stats="deletion_query_stats.csv",
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
        summary="deletion_summary.yaml",
        query_stats="deletion_query_stats.csv",
    params:
        min_deletion_size=50,
        sample_name="filtered_sample"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/deletion_block_detection"
```

### With Block Merging

Use `merge_distance` to collapse deletion blocks whose start and end positions are within a given number of bases of each other. This is useful for grouping sequencing-noise variants of the same underlying deletion into a single representative entry.

```python
rule detect_deletion_blocks:
    input:
        bams="sorted.bam"
    output:
        reads="deletion_reads.csv",
        deletions="deletion_blocks.csv",
        summary="deletion_summary.yaml",
        query_stats="deletion_query_stats.csv",
    params:
        min_deletion_size=50,
        merge_distance=10,
        sample_name="patient1"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/deletion_block_detection"
```

### Query regions

Optional `params.query` requests per-region counts in `query_stats`. Pass either:

- A single string: `"HXB2F:500-700"`
- Several regions in one string, separated by semicolons: `"HXB2F:500-700;HXB2F:800-900"`
- A list of such strings from Snakemake

Format matches `cigarmath/slice`: `ref:start-end`. Reference names may include letters, digits, `_`, `.`, and `-` (e.g. `NC_045512.2:100-200`). Coordinates use the **same overlap rule as slice** (read overlaps the interval if it would be counted as overlapping that window in the slice wrapper). Each read is attributed to `segments[0].reference_name` so queries are scoped to the correct contig in multi-reference BAMs.

**Definitions:**

- **reads_covering**: reads whose alignment span overlaps the query interval on that reference.
- **reads_with_deletion_overlapping**: reads counted above that also carry at least one reported deletion block (after `min_deletion_size`, allowedlist, and optional `merge_distance`) whose reference interval overlaps the query (half-open intervals, consistent with `cigarmath.reference_deletion_blocks`).

If `query` is absent or empty, `query_stats` is written with a header row only.

## Parameters

- `min_deletion_size` (int, optional): Minimum size of deletions to detect. Defaults to 50.
- `merge_distance` (int, optional): Maximum coordinate distance (in bases) between two deletion blocks' starts **and** ends for them to be merged into a single representative block. The representative is chosen as the block with the highest read count. Set to 0 (default) to disable merging.
- `sample_name` (str, optional): Sample name to include in metrics. Defaults to "sample".
- `query` (str, list, optional): Region(s) for `query_stats`; see [Query regions](#query-regions).

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

When `merge_distance > 0`, deletion coordinates in this file are remapped to their representative (merged) block coordinates.

### Deletion-centered CSV (`deletions`)

Contains one row per unique deletion (or merged representative when `merge_distance > 0`):

| Column | Description |
|--------|-------------|
| `deletion_start` | Start position of deletion on reference |
| `deletion_end` | End position of deletion on reference |
| `deletion_size` | Size of deletion (end - start) |
| `read_count` | Number of reads containing this deletion |
| `coverage_count` | Number of reads fully covering this deletion region |

### Query-region CSV (`query_stats`)

One row per requested region when `params.query` is non-empty; otherwise header only.

| Column | Description |
|--------|-------------|
| `region` | Original region token from `query` |
| `reference` | Reference sequence name |
| `start` | Interval start (same convention as slice) |
| `end` | Interval end |
| `reads_covering` | Reads overlapping that interval on that reference |
| `reads_with_deletion_overlapping` | Of those reads, how many have a deletion overlapping the interval |

### Summary YAML (`summary`)

Contains summary statistics for MultiQC integration:

- `sample_name`: Sample identifier
- `total_reads`: Total number of reads processed
- `reads_with_deletions`: Number of reads containing at least one deletion
- `unique_deletion_count`: Number of unique deletion blocks found (after any merging)
- `total_deletion_count`: Total number of deletions across all reads
- `deletion_frequency`: Fraction of reads with deletions
- `deletion_richness`: Number of distinct deletion blocks (equivalent to `unique_deletion_count`)
- `deletion_shannon_entropy`: Shannon entropy (nats) of the deletion frequency distribution — higher values indicate a more even spread across many deletion types
- `min_deletion_size`: Minimum deletion size threshold used
- `merge_distance`: Block-merging distance threshold used
- `input_bam_count`: Number of input BAM files processed
- `allowedlist_used`: Whether an allowedlist filter was applied
- `allowedlist_size`: Number of deletions in allowedlist (if used)

#### Diversity metrics

`deletion_richness` is the raw count of distinct deletion blocks and reflects how many different deletion events are present.

`deletion_shannon_entropy` (H) is calculated as:

```
H = -Σ p_i * ln(p_i)
```

where p_i is the proportion of all deletion observations belonging to block i. H = 0 when all reads share the same single deletion; H increases as deletions are distributed more evenly across more distinct blocks.

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
