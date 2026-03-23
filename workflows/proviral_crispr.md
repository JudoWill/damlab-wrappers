# Proviral CRISPR Pipeline

A Snakemake pipeline that automates CRISPResso2 editing analysis across one or
more samples.  It accepts reads from paired or single-end FASTQ files, or
directly from a BAM file (with optional region-level read slicing), resolves
the amplicon from either an inline sequence string or a FASTA file, and
optionally runs pairwise CRISPRessoCompare across labelled experiment/control
groups.

## Table of Contents

- [Overview](#overview)
- [Pipeline Stages](#pipeline-stages)
- [Quick Start](#quick-start)
- [samples.csv Reference](#samplescsv-reference)
- [Configuration](#configuration)
- [Input Modes](#input-modes)
- [Amplicon Resolution](#amplicon-resolution)
- [Pairwise Comparison](#pairwise-comparison)
- [Output Files](#output-files)
- [Running on a Cluster](#running-on-a-cluster)
- [Usage Examples](#usage-examples)
- [Troubleshooting](#troubleshooting)

---

## Overview

```
samples.csv
    ↓
┌──────────────────────────────────────────────┐
│  Input resolution (per sample)               │
│                                              │
│  FASTQ R1 [+ R2]  ──────────────────────┐   │
│  BAM ──→ bam2fastx ──────────────────── ┤   │
│  BAM + region ──→ cigarmath/slice ────── ┘   │
└──────────────────────────────────────────────┘
    ↓
BAM rows only: deletion_block_detection  →  deletion_detection/*  (parallel to FASTQ prep)
    ↓
CRISPResso  →  crispresso/CRISPResso_on_{sample_name}/
    ├──→  (if comparison column present)
    │     CRISPRessoCompare  →  crispresso/CRISPRessoCompare_{exp}_vs_{ctrl}/
    └──→  CRISPRessoAggregate  →  crispresso/CRISPRessoAggregate_on_all/
```

**Key features:**
- Three flexible read input modes: paired FASTQ, single-end FASTQ, or BAM
- BAM region slicing with `cigarmath/slice` — extracts only the amplicon-covering
  portion of long reads before passing them to CRISPResso
- Amplicon supplied as an inline DNA string or a path to a FASTA file
- Automatic CRISPRessoAggregate report combining all samples into one summary
- Optional automatic pairwise comparison of every experiment sample against
  every control sample using CRISPRessoCompare
- For **BAM** samples, `cigarmath/deletion_block_detection` runs on the same
  `bam_file` (read-level and block-level deletion CSVs, summary YAML, and
  optional per-region stats from `deletion_query`)
- Wrappers fetched from GitHub by default; override with a local path for
  development or reproducibility pinning

---

## Pipeline Stages

### 1. Read Preparation

Converts input reads into FASTQ for CRISPResso.  The rule used depends on the
`samples.csv` columns present for each sample (see [Input Modes](#input-modes)).

| Rule | Wrapper | When used |
|------|---------|-----------|
| `slice_bam_region` | `cigarmath/slice` | `bam_file` + `region` provided |
| `bam_to_fastq` | `cigarmath/bam2fastx` | `bam_file` only (no `region`) |
| *(none)* | — | `fastq_r1` provided directly |

Intermediate files are written to `fastq/` and are not final pipeline outputs.

### 2. Deletion block detection (`rule deletion_block_detection`)

For every sample with `bam_file` set, runs
[`cigarmath/deletion_block_detection`](../cigarmath/deletion_block_detection/)
on that BAM (independent of `slice_bam_region` / `bam_to_fastq`).

- Optional `deletion_query` column: regions in `ref:start-end` form; multiple
  regions in one cell separated by semicolons (e.g. `HXB2F:500-700;HXB2F:800-900`).
  Passed as the wrapper’s `query` parameter; see the wrapper README for column
  definitions in `deletion_query_stats.csv`.
- If `deletion_query` is blank, `deletion_query_stats.csv` still appears with a
  header-only table.
- Thresholds: `MIN_DELETION_SIZE` and `DELETION_MERGE_DISTANCE` in `run.meta.yaml`
  (defaults 50 and 10, same as the standalone deletion workflow).

Outputs under `deletion_detection/` (see [Output Files](#output-files)).

### 3. CRISPResso (`rule crispresso`)

Runs [CRISPResso2](https://github.com/pinellolab/CRISPResso2) on each sample
using the `CRISPR/crispresso-core` wrapper.

- 4 threads per sample by default
- Amplicon sequence resolved at runtime (string or FASTA — see
  [Amplicon Resolution](#amplicon-resolution))
- Guide sequence (`grna` column) passed as `--guide_seq`
- Output: `crispresso/CRISPResso_on_{sample_name}/`

### 4. CRISPRessoAggregate (`rule crispresso_aggregate`)

Runs [CRISPRessoAggregate](https://docs.crispresso.com/latest/aggregate/tool.html)
across **all** samples in the run using the `CRISPR/crispresso-aggregate` wrapper.
This step always executes regardless of whether the `comparison` column is
present.

- Input: every `CRISPResso_on_{sample_name}` directory produced in stage 3
- Output: `crispresso/CRISPRessoAggregate_on_all/` — a single HTML report and
  summary plots covering the entire run
- 4 threads

### 5. CRISPRessoCompare (`rule crispresso_compare`)

Only generated when the `comparison` column is present in `samples.csv`.
Runs `CRISPRessoCompare` for every combination of experiment × control sample
(Cartesian product).

- Input: the two `CRISPResso_on_*` directories produced in stage 3
- Output: `crispresso/CRISPRessoCompare_{exp_name}_vs_{ctrl_name}/`

---

## Quick Start

### 1. Prepare your run directory

```bash
mkdir my_crispr_run
cd my_crispr_run
```

### 2. Create `samples.csv`

**FASTQ input, amplicon as a sequence string:**

```csv
sample_name,grna,amplicon,fastq_r1
treated,TGCAGGTCGACAGATCCCCG,GCAGTCCGAAGGCTTAGATCCTGCAGGTCGACAGATCCCCGGGTACCGAG,reads/treated_R1.fastq.gz
control,TGCAGGTCGACAGATCCCCG,GCAGTCCGAAGGCTTAGATCCTGCAGGTCGACAGATCCCCGGGTACCGAG,reads/control_R1.fastq.gz
```

**BAM input with region slicing and comparison labels:**

```csv
sample_name,grna,amplicon,bam_file,region,comparison
treated_A,TGCAGGTCGACAGATCCCCG,amplicons/hxb2_target.fasta,bams/treated_A.bam,chr1:1000-1200,experiment
treated_B,TGCAGGTCGACAGATCCCCG,amplicons/hxb2_target.fasta,bams/treated_B.bam,chr1:1000-1200,experiment
untreated,TGCAGGTCGACAGATCCCCG,amplicons/hxb2_target.fasta,bams/untreated.bam,chr1:1000-1200,control
```

### 3. Create `run.meta.yaml` (optional)

If no config file is present the pipeline uses all defaults. Override settings
as needed:

```yaml
samples_csv: samples.csv
# damlab_prefix: /path/to/local/damlab-wrappers  # uncomment to use local wrappers
```

### 4. Run the pipeline

**Locally:**
```bash
cd my_crispr_run
snakemake -s /path/to/damlab-wrappers/workflows/proviral_crispr.smk \
          --use-conda --cores 8
```

**Via `data_scripts` makefile:**
```bash
# from the data_scripts directory
make proviral-crispr ROOT=/path/to/my_crispr_run MACHINE=Picotte
```

---

## samples.csv Reference

Every row is one CRISPResso run.  Either `fastq_r1` or `bam_file` must be
present.

| Column | Required | Description |
|--------|----------|-------------|
| `sample_name` | yes | Unique name for this sample. Used in all output directory names. |
| `grna` | yes | Guide RNA sequence (without PAM). Passed to CRISPResso `--guide_seq`. |
| `amplicon` | yes | Amplicon sequence string **or** path to a FASTA file containing one or more amplicon sequences. |
| `fastq_r1` | cond. | Path to R1 FASTQ (or the only FASTQ for single-end). Required unless `bam_file` is set. |
| `fastq_r2` | no | Path to R2 FASTQ for paired-end experiments. Omit or leave blank for single-end. |
| `bam_file` | cond. | Path to a BAM file. Required unless `fastq_r1` is set. |
| `region` | no | Genomic region in `chr:start-stop` format. Only used with `bam_file`. When set, reads are sliced to this region before CRISPResso. |
| `comparison` | no | `experiment` or `control`. When present, enables automatic CRISPRessoCompare runs for every experiment × control pair. |
| `deletion_query` | no | **BAM samples only.** Optional regions for per-window read/deletion counts (`ref:start-end`; multiple separated by `;`). See [Pipeline Stages](#pipeline-stages) §2. |
| `guide_name` | no | Display name for the guide RNA in CRISPResso plots (`--guide_name`). |
| `amplicon_name` | no | Display name for the amplicon (`--amplicon_name`). Defaults to the FASTA record ID when `amplicon` is a file. |
| `quantification_window_center` | no | Cleavage offset from the 3′ end of the guide sequence (`--quantification_window_center`). CRISPResso default: −3. |
| `quantification_window_size` | no | Number of bp around the cleavage site to include in quantification (`--quantification_window_size`). CRISPResso default: 1. |
| `expected_hdr_amplicon_seq` | no | Expected HDR (homology-directed repair) amplicon sequence (`--expected_hdr_amplicon_seq`). |

---

## Configuration

The pipeline reads `run.meta.yaml` from the working directory by default.
Override with `--configfile` on the command line.

| Key | Default | Description |
|-----|---------|-------------|
| `samples_csv` | `samples.csv` | Path to the samples CSV file, relative to the working directory (`ROOT`). Read with UTF-8 BOM support; column names are stripped so headers like ``deletion_query `` still match. |
| `MIN_DELETION_SIZE` | `50` | Minimum deletion length (bp) for `deletion_block_detection`. |
| `DELETION_MERGE_DISTANCE` | `10` | Merge deletion blocks whose coordinates are within this distance. |
| `DEBUG_DELETION_QUERY` | `true` | If true, the deletion wrapper writes query-param diagnostics to `logs/{sample_name}.deletion_detection.log`. Set to `false` to turn off. |
| `damlab_prefix` | GitHub `main` branch | Base URL or local path for damlab-wrappers. See below. |

### `damlab_prefix`

```yaml
# Default — fetch wrappers directly from GitHub
# (no value needed in config; this is the built-in default)

# Pin to a specific release tag
damlab_prefix: https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/tags/v1.2.3

# Use a local checkout (development / offline use)
damlab_prefix: /home/you/repos/damlab-wrappers
```

When `damlab_prefix` starts with `http://` or `https://`, the URL is used
verbatim in the `wrapper:` directive.  Otherwise it is prefixed with `file:`
so Snakemake treats it as a local path.

---

## Input Modes

Each sample is independently assigned an input mode based on which columns are
populated.

### FASTQ Mode

Provide `fastq_r1` and optionally `fastq_r2`.  Files are passed directly to
CRISPResso without any intermediate step.

```csv
sample_name,grna,amplicon,fastq_r1,fastq_r2
my_sample,ATCGATCGATCGATCGATCG,ATCG...,reads/R1.fastq.gz,reads/R2.fastq.gz
```

### BAM Mode (no region)

Provide `bam_file` without a `region`.  All mapped reads in the BAM are
exported to FASTQ using `cigarmath/bam2fastx` before being passed to
CRISPResso.  Useful when the BAM already contains only the reads of interest.

```csv
sample_name,grna,amplicon,bam_file
my_sample,ATCGATCGATCGATCGATCG,ATCG...,aligned/my_sample.bam
```

Intermediate file: `fastq/{sample_name}.bam.fastq`

### BAM + Region Mode

Provide `bam_file` and a `region` in `chr:start-stop` format.  Reads
overlapping the region are sliced by `cigarmath/slice` — only the bases
covering the target window are retained in the output FASTQ.  This is the
recommended mode when working with long-read data where each read may span far
beyond the amplicon boundaries.

```csv
sample_name,grna,amplicon,bam_file,region
my_sample,ATCGATCGATCGATCGATCG,ATCG...,aligned/my_sample.bam,chr1:2000-2250
```

Intermediate file: `fastq/{sample_name}.slice.fastq`

---

## Amplicon Resolution

The `amplicon` column accepts two forms:

**Inline sequence string** — any value that is not a path to an existing file
is treated as the amplicon DNA sequence and passed to CRISPResso via
`--amplicon_seq`.

```csv
amplicon
GCAGTCCGAAGGCTTAGATCCTGCAGGTCGACAGATCCCCGGGTACCGAGCTCGAATTC
```

**FASTA file path** — if the value is the path to an existing file it is used
as `input.amplicon_fasta`.  The FASTA record IDs become the amplicon names.
Multiple records (comma-separated in CRISPResso) are supported.

```csv
amplicon
/path/to/reference/hxb2_target_region.fasta
```

The check is performed at job execution time using `os.path.exists()`.

---

## Pairwise Comparison

When any row has a non-empty `comparison` value, the pipeline automatically
generates `CRISPRessoCompare` runs.

Rows with `comparison = experiment` and rows with `comparison = control` are
identified.  Every experiment sample is compared to every control sample
(Cartesian product).

**Example `samples.csv` with three samples generating two comparisons:**

```csv
sample_name,grna,amplicon,fastq_r1,comparison
treated_high,TGCAGGTCGACAGATCCCCG,ATCG...,reads/high.fastq.gz,experiment
treated_low,TGCAGGTCGACAGATCCCCG,ATCG...,reads/low.fastq.gz,experiment
untreated,TGCAGGTCGACAGATCCCCG,ATCG...,reads/ctrl.fastq.gz,control
```

**Comparisons generated:**
```
crispresso/CRISPRessoCompare_treated_high_vs_untreated/
crispresso/CRISPRessoCompare_treated_low_vs_untreated/
```

If the `comparison` column is absent or all values are empty/NaN, no
CRISPRessoCompare jobs are created.

---

## Output Files

```
{ROOT}/
├── samples.csv
├── run.meta.yaml
│
├── fastq/                              # Intermediate files (BAM-derived FASTQs)
│   ├── {sample_name}.bam.fastq         #   BAM mode
│   └── {sample_name}.slice.fastq       #   BAM + region mode
│
├── deletion_detection/                 # BAM samples only (cigarmath/deletion_block_detection)
│   ├── {sample_name}.deletion_reads.csv
│   ├── {sample_name}.deletion_blocks.csv
│   ├── {sample_name}.deletion_summary.yaml
│   └── {sample_name}.deletion_query_stats.csv   # header-only if no deletion_query
│
├── crispresso/
│   ├── CRISPResso_on_{sample_name}/            # CRISPResso output (one per sample)
│   │   ├── CRISPResso_output.html
│   │   ├── Alleles_frequency_table.zip
│   │   └── ...
│   ├── CRISPRessoAggregate_on_all/             # Aggregate report across all samples
│   │   ├── CRISPRessoAggregate_output.html
│   │   └── ...
│   └── CRISPRessoCompare_{exp}_vs_{ctrl}/      # Compare output (one per pair, optional)
│       ├── CRISPRessoCompare_output.html
│       └── ...
│
└── logs/
    ├── {sample_name}.crispresso.log
    ├── {sample_name}.slice.log
    ├── {sample_name}.bam2fastx.log
    ├── {sample_name}.deletion_detection.log
    ├── aggregate.log
    └── {exp_name}_vs_{ctrl_name}.compare.log
```

---

## Running on a Cluster

### Using the `data_scripts` makefile (recommended)

On shared filesystems (e.g. NFS), Snakemake may raise `MissingOutputException` right after a successful SLURM job because the submit host has not yet seen files written on a compute node. The Picotte profile sets `latency-wait: 90` (seconds). If you still see this, increase it (e.g. `snakemake ... --latency-wait 180`) or retry the same targets.

```bash
# Dry run to check the DAG
make proviral-crispr ROOT=/path/to/run MACHINE=Picotte EXTRA="-n"

# Full run on Picotte SLURM cluster
make proviral-crispr ROOT=/path/to/run MACHINE=Picotte EXTRA=""
```

The `MACHINE=Picotte` argument selects `profiles/Picotte/` which is pre-configured for the Drexel Picotte cluster:
- SLURM executor, `def` partition
- `crispresso` jobs: 4 CPUs, 16 GB RAM, 2-hour runtime
- All other rules: 1 CPU, 8 GB RAM, 4-hour default runtime

### Directly with Snakemake

```bash
snakemake \
  -s /path/to/damlab-wrappers/workflows/proviral_crispr.smk \
  -d /path/to/run \
  --configfile /path/to/data_scripts/modes/Picotte/proviral_crispr.yaml \
  --profile /path/to/data_scripts/profiles/Picotte \
  --use-conda
```

### Custom SLURM profile

```yaml
# profiles/my_cluster/config.yaml
executor: slurm
jobs: 50
use-conda: true
conda-prefix: /shared/conda

default-resources:
  slurm_account: myproject
  slurm_partition: standard
  runtime: 120
  mem_gb: 8

set-resources:
  crispresso:
    cpus_per_task: 4
    mem_gb: 16
    runtime: 120
```

---

## Usage Examples

### Example 1 — Single-end FASTQ, inline amplicon

`samples.csv`:
```csv
sample_name,grna,amplicon,fastq_r1
ctrl,TGCAGGTCGACAGATCCCCG,GCAGTCCGAAGGCTTAGATCCTGCAGGTCGACAGATCCCCGGGTACCGAGCTCGAATTC,reads/ctrl.fastq.gz
```

```bash
snakemake -s workflows/proviral_crispr.smk --use-conda --cores 4
```

Output: `crispresso/CRISPResso_on_ctrl/`

---

### Example 2 — Paired-end FASTQ, amplicon from FASTA

`samples.csv`:
```csv
sample_name,grna,amplicon,fastq_r1,fastq_r2
sample_A,TGCAGGTCGACAGATCCCCG,refs/amplicon.fasta,reads/A_R1.fastq.gz,reads/A_R2.fastq.gz
sample_B,TGCAGGTCGACAGATCCCCG,refs/amplicon.fasta,reads/B_R1.fastq.gz,reads/B_R2.fastq.gz
```

Output: `crispresso/CRISPResso_on_sample_A/`, `crispresso/CRISPResso_on_sample_B/`

---

### Example 3 — Long-read BAM with region slicing and comparison

`samples.csv`:
```csv
sample_name,grna,amplicon,bam_file,region,comparison
treated_1,TGCAGGTCGACAGATCCCCG,refs/amplicon.fasta,bams/treated_1.bam,HIV1:2550-2810,experiment
treated_2,TGCAGGTCGACAGATCCCCG,refs/amplicon.fasta,bams/treated_2.bam,HIV1:2550-2810,experiment
mock,TGCAGGTCGACAGATCCCCG,refs/amplicon.fasta,bams/mock.bam,HIV1:2550-2810,control
```

Outputs:
```
crispresso/CRISPResso_on_treated_1/
crispresso/CRISPResso_on_treated_2/
crispresso/CRISPResso_on_mock/
crispresso/CRISPRessoAggregate_on_all/
crispresso/CRISPRessoCompare_treated_1_vs_mock/
crispresso/CRISPRessoCompare_treated_2_vs_mock/
```

---

### Example 4 — Mixed input modes in one run

```csv
sample_name,grna,amplicon,fastq_r1,bam_file,region,comparison
illumina_treated,TGCAGGTCGACAGATCCCCG,ATCG...,reads/illumina.fastq.gz,,,experiment
nanopore_ctrl,TGCAGGTCGACAGATCCCCG,ATCG...,,bams/nano.bam,chr3:5000-5300,control
```

Each sample is handled independently; input mode is detected per-row.

---

## Troubleshooting

### "Sample not found in samples CSV"

The wildcard `{sample_name}` in an output path does not match any
`sample_name` value in `samples.csv`. Check for trailing whitespace or
inconsistent capitalisation in the CSV.

### Deletion `query_stats` is header-only; log shows `deletion_query` param is `None`

Snakemake is calling the wrapper with an empty query because the value read from `samples.csv` for that `sample_name` is missing or blank. The pipeline now normalizes headers (UTF-8 BOM, spaces) and matches `deletion_query` case-insensitively; `sample_name` cells are stripped.

If it still happens, confirm the **same** CSV the run uses actually contains text in `deletion_query` for that row (cluster copy vs laptop):

```bash
python -c "import pandas as pd; df=pd.read_csv('samples.csv',encoding='utf-8-sig'); print(df.columns.tolist()); print(df[['sample_name','deletion_query']])"
```

### CRISPResso fails with "no reads mapped to amplicon"

- Verify the amplicon sequence matches the reference the BAM was aligned to.
- If using `region`, confirm the coordinates are in `chr:start-stop` format and
  that the BAM contains reads mapping to that region:
  ```bash
  samtools view -c sample.bam chr1:2000-2250
  ```
- For paired-end FASTQ, check that R1 and R2 files are in the correct order.

### "Either input.amplicon_fasta or params.amplicon_seq must be provided"

The `amplicon` value in `samples.csv` was treated as a sequence string but
`amplicon_seq` ended up empty.  This can happen if the cell contains only
whitespace.  Check the CSV for blank or whitespace-only amplicon values.

### "No such file or directory" for amplicon FASTA

The `amplicon` value is a path that does not exist at the time the job runs.
The path is resolved relative to the Snakemake working directory (`ROOT`).
Use an absolute path or a path relative to `ROOT`.

### Wrapper download errors (GitHub)

By default the pipeline fetches wrappers from GitHub at each run.  If the
cluster has restricted outbound internet access, set `damlab_prefix` to a
local checkout in `run.meta.yaml`:

```yaml
damlab_prefix: /path/to/local/damlab-wrappers
```

### Rerunning failed samples

Snakemake's standard `--rerun-incomplete` flag will pick up any samples whose
output directory is missing or incomplete:

```bash
snakemake -s workflows/proviral_crispr.smk --use-conda --cores 8 --rerun-incomplete
```
