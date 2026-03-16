# BAM Split

A wrapper for dividing a SAM/BAM file into multiple output files based on one of several criteria: BAM flag attributes, auxiliary tag values, BED-file region overlap, query sequence size, or a CSV query-name→category mapping.

## Version

Current version: 1.0.0

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

### Split by Flag

```python
rule split_mapped_unmapped:
    input:
        bam="aligned.bam"
    output:
        directory("split_flag/")
    params:
        split_by="flag",
        flag_categories={"mapped": {"is_unmapped": False}, "unmapped": {"is_unmapped": True}},
        output_format="fasta"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/split"
```

### Split by Auxiliary Tag

```python
rule split_by_haplotype:
    input:
        bam="phased.bam"
    output:
        dir=directory("split_hp/"),
        summary="split_summary.yaml"
    params:
        split_by="tag",
        tag_name="HP",
        output_format="bam",
        top_n=2
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/split"
```

### Split by BED Regions

The 4th column of the BED file is used as the category name.

```python
rule split_by_region:
    input:
        bam="aligned.bam",
        bed_file="regions.bed"
    output:
        directory("split_regions/")
    params:
        split_by="bed",
        output_format="fastq"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/split"
```

### Split by Query Size

```python
rule split_by_size:
    input:
        bam="aligned.bam"
    output:
        directory("split_size/")
    params:
        split_by="query_size",
        size_bins={"short": [0, 500], "medium": [500, 5000], "long": [5000, null]},
        output_format="fastq"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/split"
```

### Split by CSV Query-Name Mapping

The CSV must have at least two columns; the first is `query_name` and the second is the category.

```python
rule split_by_label:
    input:
        bam="aligned.bam",
        csv_file="read_labels.csv"
    output:
        directory("split_csv/")
    params:
        split_by="csv",
        output_format="bam"
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/split"
```

### With Required Reference Overlap Filter

Only reads that fully span **all** listed regions (as a single uninterrupted mapping block) are passed to the split criterion. Reads that fail this filter are counted but not written to any output file.

```python
rule split_nfl_haplotype:
    input:
        bam="aligned.bam"
    output:
        dir=directory("split/"),
        summary="split_summary.yaml"
    params:
        split_by="tag",
        tag_name="HP",
        output_format="bam",
        required_reference_overlap=["HIV:790-2292", "HIV:4920-6225"]
    wrapper:
        "file:path/to/damlab-wrappers/cigarmath/split"
```

## Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `split_by` | str | **required** | One of `"flag"`, `"tag"`, `"bed"`, `"query_size"`, `"csv"` |
| `output_format` | str | **required** | One of `"fasta"`, `"fastq"`, `"sam"`, `"bam"` |
| `flag_categories` | dict | — | `{category_name: {pysam_attr: bool, ...}}` — first matching category wins. Attributes are pysam `AlignedSegment` boolean properties (e.g. `is_unmapped`, `is_supplementary`, `is_secondary`, `is_duplicate`). |
| `tag_name` | str | — | BAM auxiliary tag to split on (e.g. `"HP"`, `"RG"`). |
| `tag_categories` | dict | — | Optional `{raw_tag_value: category_name}` map. If omitted, the raw string value of the tag is used as the category name. |
| `size_bins` | dict | — | `{category_name: [min, max]}` — first matching bin wins. Use `null` for an open-ended bound. |
| `top_n` | int | `null` | When set, performs a pre-pass to count categories and keeps only the top-N by read count. Reads in lower-ranked categories are routed to `unclassified`. |
| `required_reference_overlap` | list[str] | `null` | List of `"chr:start-stop"` region strings (0-based, half-open). A read must have a single uninterrupted mapping block that fully spans **every** listed region. Reads that fail are silently excluded (not written to any output file, including `unclassified`). |
| `sample_name` | str | `"sample"` | Sample identifier written to the summary YAML. |
| `unclassified_name` | str | `"unclassified"` | File stem used for reads that match no category. |

## Input

- `bam` (positional): SAM or BAM file. The wrapper reads the `SO` header field to detect name-sorted/collated BAMs and uses `cm.io.combined_segment_stream` for those, merging supplementary alignments before applying any filter or split criterion.
- `bed_file` (named, optional): BED file required when `split_by="bed"`.
- `csv_file` (named, optional): CSV file required when `split_by="csv"`.

## Output

- `output[0]` (positional): Output directory. The wrapper creates one file per category using the naming pattern `{category}.{ext}`, where `ext` depends on `output_format`.
- `output[1]` (positional, optional): Summary YAML path.

### Name-Sorted / Collated BAMs and Supplementary Alignments

When the input BAM header declares `SO:queryname` or `SO:collated`, the wrapper groups all alignment records sharing the same `query_name` (primary + supplementary) into a single logical unit via `cm.io.combined_segment_stream`. All filter and split logic is then applied to the merged combined CIGAR:

- **`required_reference_overlap`**: evaluated against the merged CIGAR, so a chimeric read whose two halves together span a required region will pass even if neither half alone spans it.
- **`query_size`**: computed from the merged query block, reflecting the full read length.
- **SAM/BAM output**: all segments in the group are written to the output file, preserving supplementary records.
- **FASTA/FASTQ output**: only `segments[0].query_sequence` is written (one record per read group).

For coordinate-sorted BAMs each record is evaluated individually.

### SAM/BAM Header

When `output_format` is `"sam"` or `"bam"`, each output file receives a copy of the input header with an additional `PG` entry:

```
@PG  ID:cigarmath-split  PN:cigarmath-split  VN:1.0.0  PP:<previous_PG_ID>
```

### Summary YAML

```yaml
# Cigarmath Split
sample_name: patient1
split_by: tag
tag_name: HP
output_format: bam
output_directory: results/split/
is_name_sorted: true
top_n: null
required_reference_overlap:
  - "HIV:790-2292"
  - "HIV:4920-6225"
total_reads: 5000
required_overlap_filtered: 312
categories:
  HP_0: 2500
  HP_1: 2400
unclassified: 100
```

## Flag Attribute Reference

Common pysam `AlignedSegment` boolean attributes usable in `flag_categories`:

| Attribute | Description |
|---|---|
| `is_unmapped` | Read is unmapped |
| `is_mapped` | Read is mapped (alias) |
| `is_paired` | Read is paired-end |
| `is_proper_pair` | Both mates mapped in expected orientation |
| `is_read1` | First mate in pair |
| `is_read2` | Second mate in pair |
| `is_secondary` | Secondary alignment |
| `is_supplementary` | Supplementary (chimeric) alignment |
| `is_duplicate` | PCR or optical duplicate |
| `is_qcfail` | Fails platform/vendor QC |

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors

* Will Dampier, PhD

## Software Requirements

* [pysam](https://pysam.readthedocs.io/) (tested with v0.19.0)
* [cigarmath](https://github.com/DamLabResources/cigarmath) (tested with v0.1.0)
