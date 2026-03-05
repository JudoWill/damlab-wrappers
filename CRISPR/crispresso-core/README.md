# CRISPResso Wrapper

Version: 1.0.0

Wraps [CRISPResso](https://docs.crispresso.com/latest/core/tool.html) (v2.3+) for CRISPR editing quantification on a single amplicon.

## Input

* `fastq_r1` (required): FASTQ file of reads (R1 for paired-end, or single-end reads).
* `fastq_r2` (optional): FASTQ file of R2 reads for paired-end experiments.
* `amplicon_fasta` (optional): FASTA file containing one or more amplicon sequences. Mutually exclusive with `params.amplicon_seq`. Record IDs are used as amplicon names unless `params.amplicon_name` is also set.

## Output

* `directory("results/CRISPResso_on_{name}")`: the CRISPResso output directory. The directory name **must** follow the `CRISPResso_on_{name}` convention, or `params.name` must be set explicitly so the wrapper can derive the correct `--name` and `--output_folder` arguments to pass to CRISPResso.

## Parameters

* `amplicon_seq` (str, required if `amplicon_fasta` not given): Amplicon sequence(s), comma-separated for multiple amplicons.
* `amplicon_name` (str, optional): Amplicon name(s), comma-separated. Defaults to FASTA record IDs when `amplicon_fasta` is used.
* `name` (str, optional): Run name. Defaults to the directory basename with `CRISPResso_on_` stripped.
* `guide_seq` (str, optional): sgRNA sequence(s), comma-separated.
* `guide_name` (str, optional): sgRNA name(s), comma-separated.
* `expected_hdr_amplicon_seq` (str, optional): Expected HDR amplicon sequence.
* `quantification_window_size` (int, optional): Number of bp around the cleavage site to quantify. CRISPResso default: 1.
* `quantification_window_center` (int, optional): Cleavage offset from 3' end of guide. CRISPResso default: -3.
* `min_average_read_quality` (int, optional): Minimum average read quality to keep a read. CRISPResso default: 0.
* `min_single_bp_quality` (int, optional): Minimum quality at any single base to keep a read. CRISPResso default: 0.
* `ignore_substitutions` (bool, optional): Ignore substitution events. Default: False.
* `ignore_insertions` (bool, optional): Ignore insertion events. Default: False.
* `ignore_deletions` (bool, optional): Ignore deletion events. Default: False.
* `discard_indel_reads` (bool, optional): Discard reads with indels. Default: False.
* `base_editor_output` (bool, optional): Produce base editor output. Default: False.
* `conversion_nuc_from` (str, optional): Base to convert from (for base editor analysis).
* `conversion_nuc_to` (str, optional): Base to convert to (for base editor analysis).
* `extra` (str, optional): Any additional CRISPResso arguments passed verbatim.

## Examples

### Amplicon sequence as a parameter

```python
rule crispresso:
    input:
        fastq_r1 = "data/sample_R1.fastq.gz",
        fastq_r2 = "data/sample_R2.fastq.gz",
    output:
        directory("results/CRISPResso_on_sample")
    params:
        amplicon_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        guide_seq    = "ATCGATCGATCGATCGATCG",
    threads: 4
    log: "logs/crispresso_sample.log"
    wrapper:
        "file:path/to/damlab-wrappers/CRISPR/crispresso-core"
```

### Amplicon sequence from a FASTA file

```python
rule crispresso:
    input:
        fastq_r1      = "data/sample_R1.fastq.gz",
        amplicon_fasta = "references/amplicon.fasta",
    output:
        directory("results/CRISPResso_on_sample")
    params:
        name      = "sample",
        guide_seq = "ATCGATCGATCGATCGATCG",
    threads: 4
    log: "logs/crispresso_sample.log"
    wrapper:
        "file:path/to/damlab-wrappers/CRISPR/crispresso-core"
```

## Software Requirements

* [CRISPResso2](https://github.com/pinellolab/CRISPResso2) ≥ 2.3 (installed via bioconda: `crispresso2`)

## Author

Will Dampier — wnd22@drexel.edu
