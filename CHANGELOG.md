# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## Unreleased

Collect wrapper level changes here until merged.

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