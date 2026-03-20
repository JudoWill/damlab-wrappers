# samtools/consensus

Wrapper for [`samtools consensus`](https://www.htslib.org/doc/samtools-consensus.html).

Generates a single FASTA consensus sequence from a coordinate-sorted BAM file.

## Inputs

| Name | Description |
|------|-------------|
| (positional 0) | Coordinate-sorted BAM file |

## Outputs

| Name | Description |
|------|-------------|
| (positional 0) | FASTA file containing the consensus sequence |

## Parameters

| Name | Default | Description |
|------|---------|-------------|
| `mode` | `"bayesian"` | Consensus mode: `"bayesian"` or `"simple"` |
| `extra` | `""` | Extra flags passed verbatim to `samtools consensus` |
| `max_reads` | omitted / `None` / `0` | When set to a positive integer `M`, alignments are counted with `samtools view -c -F <count_filter_flags>`. If the count exceeds `M`, the BAM is subsampled with `samtools view -b --subsample <M/N> --subsample-seed <subsample_seed>` into a file under the job tmpdir, then consensus runs on that file. If the count is ≤ `M`, consensus runs on the original BAM (no extra `view` pass). |
| `subsample_seed` | `0` | Passed to `samtools view --subsample-seed` when downsampling. |
| `count_filter_flags` | `"0x900"` | Exclude mask for the count step (`samtools view -c -F ...`). Default excludes secondary (`0x100`) and supplementary (`0x800`) alignments. |

## Example rule

```python
rule consensus:
    input:
        "aligned/{sample}.sorted.bam"
    output:
        "results/{sample}.consensus.fa"
    params:
        mode="bayesian",
        max_reads=50000,  # optional cap before consensus
        subsample_seed=0,
        count_filter_flags="0x900",
    threads: 4
    log:
        "logs/{sample}.consensus.log"
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/samtools/consensus/"
```

## Example rule (no downsampling)

Omit `max_reads` (or set `0`) for the original single-step `samtools consensus` behavior.

```python
rule consensus:
    input:
        "aligned/{sample}.sorted.bam"
    output:
        "results/{sample}.consensus.fa"
    params:
        mode="bayesian",
    threads: 4
    log:
        "logs/{sample}.consensus.log"
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/samtools/consensus/"
```
