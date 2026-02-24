# BAM to FASTA/FASTQ Wrapper

This wrapper converts BAM/SAM files to FASTA or FASTQ format using the cigarmath library. It supports filtering by mapped reads only and by genomic region.

## Input

- BAM or SAM file

## Output

- FASTA file (`.fa`/`.fasta`) or FASTQ file (`.fq`/`.fastq`), inferred from the output path extension

## Parameters

- `mapped_only` (bool, optional): If `True`, only output reads that are mapped. Default: `False`.
- `region` (str, optional): Genomic region to extract (e.g. `"chr1:1000-2000"` or `"chr1"`). For BAM files, requires an index (`.bai`). For SAM, region filtering is not supported.
- `output_format` (str, optional): Override output format: `"fasta"` or `"fastq"`. By default, format is inferred from the output file extension.
- `min_mapq` (int, optional): Minimum mapping quality. Default: `0`.

## Example Usage

```python
rule bam2fastx:
    input: "aligned.bam"
    output: "reads.fastq"
    params:
        mapped_only=True,
        region="chr1:1000-2000"
    wrapper: "damlab-wrappers/cigarmath/bam2fastx"
```

## Output Format

- **FASTA**: Sequence wrapped at 60 characters per line.
- **FASTQ**: Uses `query_qualities` when present. If absent (e.g. some SAM files), uses placeholder quality `!` (Phred 0) so the output is valid FASTQ.

## Notes

- Region filtering with BAM requires an index file (`.bai`). Without an index, the wrapper iterates over all reads.
- Uses cigarmath for streaming; no samtools CLI dependency.
