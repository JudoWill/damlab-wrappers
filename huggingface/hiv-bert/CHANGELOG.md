# Changelog

All notable changes to the HIV-BERT wrapper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.2] - 2025-05-15
  - Changed testing strategy for new test harness.


## [1.0.1] - 2025-04-23
 - Fixed `environment.yaml` to properly include pillow
 - Adjusted `env.yaml` to be a minimal set
 - Added `gpu` resources to the test rules to ensure it plays nice with multiple cores.


## [1.0.0] - 2025-04-10

### Added
- Initial release of the HIV-BERT wrapper
- Support for multiple HIV-BERT models:
  - damlab/hiv_bert (embedding model)
  - damlab/HIV_V3_bodysite (classification model)
  - damlab/HIV_V3_coreceptor (classification model)
- Support for multiple input formats:
  - FASTA files
  - FASTQ files
  - BAM/SAM files
- Automatic DNA to amino acid translation
- Model caching functionality
- GPU support with automatic detection
- Comprehensive error handling
- Detailed logging
- Support for custom environments
- Batch processing of sequences
- Sequence length filtering
- Output metrics generation
- Support for both embedding and classification models
- Automatic handling of rare amino acids
- Support for multiple reading frames
- BAM/SAM specific features (mapped reads only option)

### Changed
- None

### Deprecated
- None

### Removed
- None

### Fixed
- None

### Security
- None 