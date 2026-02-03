# Changelog

All notable changes to this wrapper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.1] - 2025-05-14
  - Changed testing strategy for new test harness.

## [1.0.0] - 2025-04-10

### Added
- First stable release of the tree rerooting wrapper
- Basic tree rerooting functionality using DendroPy
- Support for tree manipulation parameters:
  - schema: Tree file format
  - preserve_branch_lengths: Branch length preservation
  - preserve_support_values: Support value preservation
- Version checking functionality
- Input file validation
- Comprehensive error handling
- Detailed documentation in README.md
- Type hints for better code maintainability

### Changed
- Improved error handling with specific error messages
- Better organized code structure
- Enhanced parameter handling with sensible defaults
- Added preservation of underscores in taxon names
- Added bipartition updating during rerooting

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A

## [0.0.0] - 2025-04-10

### Added
- Initial development version
- Basic structure and functionality
- Preliminary implementation of tree rerooting wrapper 