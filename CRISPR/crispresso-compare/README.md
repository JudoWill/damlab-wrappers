# CRISPRessoCompare Wrapper

Version: 1.0.0

Wraps [CRISPRessoCompare](https://docs.crispresso.com/latest/compare/tool.html) (v2.3+) to compare editing outcomes between two CRISPResso runs — for example, a treated sample versus an untreated control.

## Input

* `folder_1` (required): Path to the first CRISPResso output directory (e.g. treated sample).
* `folder_2` (required): Path to the second CRISPResso output directory (e.g. control sample). Both runs must have used the same reference amplicon and settings.

## Output

* `directory("results/CRISPRessoCompare_{name}")`: the CRISPRessoCompare output directory. The directory name **must** follow the `CRISPRessoCompare_{name}` convention, or `params.name` must be set explicitly so the wrapper can derive the correct `--name` and `--output_folder` arguments.

## Parameters

* `name` (str, optional): Run name. Defaults to the directory basename with `CRISPRessoCompare_` stripped.
* `sample_1_name` (str, optional): Display name for sample 1 in plots.
* `sample_2_name` (str, optional): Display name for sample 2 in plots.
* `reported_qvalue_cutoff` (float, optional): Q-value cutoff for reporting significant differences. CRISPResso default: 0.05.
* `min_frequency_alleles_around_cut_to_plot` (float, optional): Minimum allele frequency to include in plots. CRISPResso default: 0.2.
* `max_rows_alleles_around_cut_to_plot` (int, optional): Maximum number of allele rows to display in plots.
* `suppress_report` (bool, optional): Suppress HTML report generation. Default: False.
* `extra` (str, optional): Any additional CRISPRessoCompare arguments passed verbatim.

## Example

```python
rule crispresso_core_treated:
    input:
        fastq_r1 = "data/treated_R1.fastq.gz",
    output:
        directory("results/CRISPResso_on_treated"),
    params:
        amplicon_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        guide_seq    = "ATCGATCGATCGATCGATCG",
    threads: 4
    log: "logs/crispresso_treated.log"
    wrapper:
        "file:path/to/damlab-wrappers/CRISPR/crispresso-core"


rule crispresso_core_control:
    input:
        fastq_r1 = "data/control_R1.fastq.gz",
    output:
        directory("results/CRISPResso_on_control"),
    params:
        amplicon_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        guide_seq    = "ATCGATCGATCGATCGATCG",
    threads: 4
    log: "logs/crispresso_control.log"
    wrapper:
        "file:path/to/damlab-wrappers/CRISPR/crispresso-core"


rule crispresso_compare:
    input:
        folder_1 = "results/CRISPResso_on_treated",
        folder_2 = "results/CRISPResso_on_control",
    output:
        directory("results/CRISPRessoCompare_treated_vs_control"),
    params:
        sample_1_name = "treated",
        sample_2_name = "control",
    log: "logs/crispresso_compare.log"
    wrapper:
        "file:path/to/damlab-wrappers/CRISPR/crispresso-compare"
```

## Software Requirements

* [CRISPResso2](https://github.com/pinellolab/CRISPResso2) ≥ 2.3 (installed via bioconda: `crispresso2`)

## Author

Will Dampier — wnd22@drexel.edu
