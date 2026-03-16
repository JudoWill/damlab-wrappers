# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-03-16

### Added

- Initial release of `cigarmath/split` wrapper
- Split BAM/SAM files by BAM flag attributes (`split_by="flag"`)
- Split by supplemental BAM tag value (`split_by="tag"`)
- Split by BED file region overlap (`split_by="bed"`)
- Split by query sequence size (`split_by="query_size"`)
- Split by CSV query-name→category mapping (`split_by="csv"`)
- `required_reference_overlap` parameter to pre-filter reads that do not fully span all specified regions; filtering is evaluated on the post-merge combined CIGAR for name-sorted BAMs
- `top_n` parameter for a two-pass pre-count and top-N category selection
- Automatic detection of name-sorted/collated BAMs via header `SO` field; uses `cm.io.combined_segment_stream` for merged supplementary alignment handling when sorted
- Output format selection: `fasta`, `fastq`, `sam`, `bam`
- SAM/BAM outputs include a copy of the input header with an appended `PG` line
- Optional summary YAML output with read counts, category tallies, and filter statistics
