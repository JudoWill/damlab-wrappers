# Changelog

All notable changes to the Pandas Merge Wrapper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.2] - 2025-05-15
  - Changed testing strategy for new test harness.

## [1.1.1] - 2025-04-30

### Fixed
 - Fixed a bug where merge-keys keep getting renamed between merges when merging with multiple keys across multiple samples.


## [1.1.0] - 2025-04-10

### Added
- Ability to merge tables with different column names for their keys
- Better error handling and logging for merge operations

### Changed
- Improved documentation with examples for different column name merging
- Updated type hints to better reflect the new functionality

## [1.0.0] - 2025-04-09

### Added
- Initial release with basic merge functionality
- Support for multiple file merging with suffixes
- Support for different merge types (inner, outer, left, right)
- Basic documentation and examples 