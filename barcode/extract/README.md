# Barcode and UMI Extractor

A wrapper for extracting cellular barcodes and UMIs from BAM files using regex patterns. This wrapper uses the regex library to accurately identify and extract barcode and UMI sequences from sequencing reads.

## Version

Current version: 1.1.0

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../README.md).

### Using Built-in Pattern
```python
rule extract_barcodes:
    input:
        "input.bam"
    output:
        "output.bam",
        "metrics.yaml"  # Optional metrics output
    params:
        builtin="SIVmac239m2",
        barcode_tag="CR",
        umi_tag="OX",
        mapped_only=False,
        sample_name="sample1"
    wrapper:
        "file:path/to/damlab-wrappers/barcode/extract"
```

### Using Custom Patterns
```python
rule extract_barcodes:
    input:
        "input.bam"
    output:
        "output.bam",
        "metrics.yaml"  # Optional metrics output
    params:
        left_primer="(PRIMER1){i<=2}(?<umi>.{18})(PRIMER2){i<=2}",
        right_primer="(PRIMER3){i<=2}(?<umi>.{18})(PRIMER4){i<=2}",
        barcode_primer="(PRIMER5){i<=2}(?<umi>.{34})(PRIMER6){i<=2}",
        barcode_tag="CR",
        umi_tag="OX",
        mapped_only=False,
        sample_name="sample1"
    wrapper:
        "file:path/to/damlab-wrappers/barcode/extract"
```

## Parameters

- `barcode_tag` (str, optional): SAM tag to use for cellular barcode sequence. Defaults to 'CR'.
- `umi_tag` (str, optional): SAM tag to use for UMI sequence. Defaults to 'OX'.
- `mapped_only` (bool, optional): Only process mapped reads. Defaults to False.
- `sample_name` (str, optional): Sample name to include in metrics. Defaults to None.

Either specify a built-in pattern:
- `builtin` (str, optional): Name of built-in primer pattern set to use. Available options:
  - `SIVNFL` : Primers for the current Nanopore SIV-NFL PCR strategy
  - Any of from [Figure 5](https://journals.asm.org/doi/full/10.1128/jvi.01420-19) :
    - `SIVmac239M`
    - `SIVmac239m2`
    - `SIVmac239Opt5M`
    - `bSHIV1054M`
    - `bSHIVAD8EOM`
    - `cSHIV174M`
    - `cSHIV224M`
  - Delimit multiple builtins with commas, ie `builtin = 'SIVNFL,SIVmac239m2'` to load both together.

Or specify custom patterns:
- `left_primer` (str, required if builtin not specified): Regex pattern for left primer with UMI capture group
- `right_primer` (str, required if builtin not specified): Regex pattern for right primer with UMI capture group
- `barcode_primer` (str, required if builtin not specified): Regex pattern for barcode primer with UMI capture group

## Input
* BAM file containing sequencing reads

## Output
* BAM file with added tags for barcodes and UMIs
* YAML metrics file (optional) containing:
  - Total reads processed
  - Number of mapped reads
  - Extraction success counts for barcodes and UMIs
  - Pairwise success counts
  - Sample name (if provided)

## Tags Added
* CR (default): Cellular barcode sequence (34bp)
* OX (default): Combined left and right UMI sequences (36bp total)

## Error Handling

The wrapper includes error handling for:
- Invalid BAM files
- Missing or invalid regex patterns
- Invalid built-in pattern names
- Missing required parameters

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors
* Will Dampier, PhD

## Software Requirements
* [pysam](https://pysam.readthedocs.io/) (tested with v0.19.0)
* [regex](https://pypi.org/project/regex/) (tested with v2023.10.3) 