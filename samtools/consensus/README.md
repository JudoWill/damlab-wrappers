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

## Example rule

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
