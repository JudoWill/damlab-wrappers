# Changelog

All notable changes to this wrapper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.1] - 2025-05-14
  - Changed testing strategy for new test harness.

## [1.1.0] - 2025-05-1

### Fixed
 - The counting algorithm conflated deletions _within_ and deletions of the _entire_ region. Adjusted outputs to split into `full_deleted` and `partial_deleted`.
 - Improved looping logic along deletion blocks. Much faster.


## [1.0.0] - 2025-04-09

### Added
- First stable release of the cigarmath/deletion_frequency wrapper
- Basic deletion frequency calculation functionality
- Support for region-specific analysis
- Support for minimum deletion size filtering
- Support for read-level statistics output
- Comprehensive metrics generation
- Detailed documentation in README.md

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

## [0.0.0] - 2024-04-09

### Added
- Initial development version
- Basic structure and functionality
- Preliminary implementation of deletion frequency calculation features 