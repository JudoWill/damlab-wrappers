# Changelog

All notable changes to this wrapper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.1] - 2025-02-24

- Initial release
- Convert BAM/SAM to FASTA or FASTQ
- Support `mapped_only` and `region` parameters
- Use cigarmath for streaming (no samtools CLI)
