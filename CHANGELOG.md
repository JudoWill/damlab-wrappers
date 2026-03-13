# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## Unreleased

Collect wrapper level changes here until merged.

### Added
  - [`hiv/vadr-genbank`](hiv/vadr-genbank/README.md) [1.0.0] - New wrapper that normalises NCBI GenBank flat files and FASTA files for use with `v-build.pl`. Rewrites the LOCUS name token to match the bare accession and strips `.N` version suffixes from the VERSION line. The FASTA is copied unchanged (VADR requires the `accession.N` header format). Includes unit tests.
  - [`hiv/vadr-build`](hiv/vadr-build/README.md) [1.0.0] - New wrapper for `v-build.pl` to build a VADR homology model from a GenBank accession. Supports online NCBI fetch (accession only) or fully offline local-file build (via `hiv/vadr-genbank` pre-processing). Includes integration tests.
  - [`hiv/vadr-annotate`](hiv/vadr-annotate/README.md) [1.0.0] - New wrapper for `v-annotate.pl`. Auto-detects `--mkey` from the model directory. `mode='hiv'` pre-configures `--alt_pass` for all recoverable HIV alert codes so that proviral sequences with frameshifts, internal stops, and other expected defects are classified as passing. Includes integration tests.
  - [`hiv/vadr-tbl2gbk`](hiv/vadr-tbl2gbk/README.md) [1.0.0] - New wrapper for `table2asn`. Accepts a `vadr-annotate` output directory and converts passing sequences + NCBI 5-column feature tables to annotated GenBank flat files (`.gbf`) and ASN.1 submission files (`.sqn`). Handles multi-sequence inputs by splitting, processing, and concatenating per-sequence. Includes integration tests.
  - `workflows/proviral_reconstruction` [new] - New end-to-end pipeline for proviral haplotype reconstruction from Nanopore data. Steps: FASTQ/BAM input → optional filtlong downsampling → Strainline haplotype reconstruction → clipqs clipping/orientation → optional VADR annotation (three modes: pre-built model, NCBI fetch, or fully offline local files) → GenBank flat file generation → optional per-cohort MSA/FastTree/phytreeviz phylogenetics. See [`workflows/proviral_reconstruction.md`](workflows/proviral_reconstruction.md).
  - [`CRISPR/crispresso-core`](CRISPR/crispresso-core/README.md) [1.0.0] - Wrapper for CRISPResso2 CRISPR editing quantification on a single amplicon. Amplicon sequence may be supplied as an inline string parameter or as a FASTA input file.
  - [`CRISPR/crispresso-compare`](CRISPR/crispresso-compare/README.md) [1.0.0] - Wrapper for CRISPRessoCompare pairwise comparison of two CRISPResso output directories.
  - [`CRISPR/crispresso-aggregate`](CRISPR/crispresso-aggregate/README.md) [1.0.0] - Wrapper for CRISPRessoAggregate multi-run aggregation. Combines any number of CRISPResso output directories into a single HTML report and summary plots.
  - [`CRISPR/crispresso-aggregate`](CRISPR/crispresso-aggregate/README.md) [1.0.0] - Wrapper for CRISPRessoAggregate to combine any number of CRISPResso runs into a single summary report.
  - [`cigarmath/bam2fastx`](cigarmath/bam2fastx/README.md) [0.0.4] - Added `primary_only` param to skip secondary and supplementary alignments. Added `unique_only` param to emit only the first occurrence of each `query_name`. Added `min_query_length` and `max_query_length` params for read-length filtering.
  - [`cigarmath/slice`](cigarmath/slice/README.md) [1.0.0] - Wrapper for extracting reads overlapping a genomic region from a BAM, slicing each read to return only bases covering the target window.
  - [`cigarmath/bam2csv`](cigarmath/bam2csv/README.md) [0.0.1] - Wrapper for extracting fields from a SAM/BAM file into a CSV.
  - `workflows/proviral_crispr` [new] - End-to-end CRISPResso2 automation pipeline. Accepts paired/single-end FASTQ or BAM input (with optional region slicing for long reads), resolves amplicons from sequence strings or FASTA files, always produces a CRISPRessoAggregate report across all samples, and optionally runs CRISPRessoCompare across labelled experiment/control groups. See [`workflows/proviral_crispr.md`](workflows/proviral_crispr.md).
  - `workflows/proviral_nfl` - Added strainline haplotype reconstruction and deletion block detection to the pipeline
  - New `rules/strainline.smk` module with `bam_to_fasta` and `strainline` rules (outputs to `strainline/` directory)
  - New `rules/deletion_detection.smk` module with `deletion_block_detection` rule (outputs to `deletion_detection/` directory)
  - Analysis outputs now included in MultiQC report
  - New config options: `STRAINLINE_PREFIX`, `MIN_DELETION_SIZE`

## [0.0.2] - 2025-05-19

### Added
 - `cigarmath/bam2gel` [1.0.1] - Added `log_scale` to  to improve resolution in gel images.
 - [`barcode/extract`](barcode/extract/README.md) [1.1.0] - Added builtin support for more SIV barcoded viruses.

## Fixed
 - `cigarmath/deletion_frequency` [1.1.0] - Split full and partial deletions into different fields for improved downstream visualization.

### Docs
  - [`cigarmath/bam2gel`](cigarmath/bam2gel/README.md) - Improved customization documentation for gels.

### Testing
  - Improved testing harness that `include`s all tests into a single Snakefile allowing all tests to be run in the same DAG. Significantly improves test speed on clusters.


## [0.0.1] - 2025-05-07

Starting gobal changelog.
