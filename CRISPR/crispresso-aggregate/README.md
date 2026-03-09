# CRISPRessoAggregate Wrapper

Version: 1.0.0

Wraps [CRISPRessoAggregate](https://docs.crispresso.com/latest/aggregate/tool.html) (v2.3+) to combine output from any number of CRISPResso runs into a single summary report.

## Input

* `crispresso_dirs` (required): One or more CRISPResso output directories to aggregate. Supply as a named list for clarity, or as positional inputs.

## Output

* `directory("results/CRISPRessoAggregate_on_{name}")`: the aggregate output directory. The directory name **must** follow the `CRISPRessoAggregate_on_{name}` convention, or `params.name` must be set explicitly.

## Parameters

* `name` (str, optional): Run name. Defaults to the directory basename with `CRISPRessoAggregate_` stripped.
* `min_reads_for_inclusion` (int, optional): Minimum number of reads a run must have to be included in the summary.
* `max_samples_per_summary_plot` (int, optional): Maximum samples per page in the PDF summary plot. Keep below ~150 to avoid memory issues.
* `suppress_report` (bool, optional): Suppress HTML report generation. Default: False.
* `suppress_plots` (bool, optional): Suppress plot generation. Default: False.
* `extra` (str, optional): Any additional CRISPRessoAggregate arguments passed verbatim.

## Example

```python
rule crispresso_aggregate:
    input:
        crispresso_dirs = expand(
            "results/CRISPResso_on_{sample}",
            sample = ["treated", "control", "mock"],
        ),
    output:
        directory("results/CRISPRessoAggregate_on_all_samples"),
    params:
        min_reads_for_inclusion = 100,
    threads: 4
    log: "logs/crispresso_aggregate.log"
    wrapper:
        "file:path/to/damlab-wrappers/CRISPR/crispresso-aggregate"
```

Each input directory is passed as a separate `--prefix` argument to CRISPRessoAggregate, so the wrapper accepts directories from any location — they do not need to share a common parent.

## Software Requirements

* [CRISPResso2](https://github.com/pinellolab/CRISPResso2) ≥ 2.3 (installed via bioconda: `crispresso2`)

## Author

Will Dampier — wnd22@drexel.edu
