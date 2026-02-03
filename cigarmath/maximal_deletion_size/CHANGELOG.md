# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.1] - 2025-05-14
  - Changed testing strategy for new test harness.

## [1.0.0] - 2025-04-28

### Added
- Initial release of the maximal deletion size wrapper
- Basic functionality to calculate reference block size and largest deletion for each read
- Output of both YAML metrics and CSV with per-read information
- Chromosome filtering capability
- Comprehensive test suite

### Features
- Processes BAM files sorted by read name
- Groups aligned segments of the same query name
- Calculates reference block size and largest deletion for each read
- Outputs summary metrics in YAML format
- Outputs detailed per-read information in CSV format
- Optional chromosome filtering
- Accurate read counting including unmapped reads 