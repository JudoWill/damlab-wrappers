# Changelog

All notable changes to this wrapper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-03-13

### Added
- `merge_distance` parameter to merge deletion blocks with similar start/end coordinates into a single representative block using union-find grouping
- `deletion_richness` field in summary YAML (count of distinct deletion blocks after merging)
- `deletion_shannon_entropy` field in summary YAML (Shannon entropy of deletion frequency distribution across blocks)
- `merge_distance` field in summary YAML reflecting the parameter used
- `calculate_shannon_entropy` helper function

### Changed
- `merge_blocks` now returns a tuple of `(Counter, dict)` where the dict maps each original block to its representative, enabling read-level deletion coordinate remapping
- Version bumped to 1.1.0

## [1.0.0] - 2026-02-23

### Added
- First stable release of the cigarmath/deletion_block_detection wrapper
- Support for multiple BAM file inputs processed together
- Read-centered CSV output with deletion information per read
- Deletion-centered CSV output with read counts and coverage statistics
- Summary YAML output for MultiQC integration
- Optional allowedlist parameter to filter specific deletions
- Configurable minimum deletion size threshold
- Comprehensive documentation in README.md

### Changed
- N/A

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A
