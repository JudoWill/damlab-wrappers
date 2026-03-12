# Proviral Reconstruction Pipeline

A Snakemake pipeline that automates long-read proviral haplotype reconstruction
from nanopore sequencing data.  It accepts reads from FASTQ or BAM files,
optionally downsamples with filtlong, reconstructs haplotypes with Strainline,
clips and orients them against a reference, and optionally runs VADR annotation
and per-cohort phylogenetic tree generation.

## Table of Contents

- [Overview](#overview)
- [Pipeline Stages](#pipeline-stages)
- [Quick Start](#quick-start)
- [samples.csv Reference](#samplescsv-reference)
- [Configuration](#configuration)
- [Annotation Modes](#annotation-modes)
- [Phylogenetic Trees](#phylogenetic-trees)
- [Output Files](#output-files)
- [Running on a Cluster](#running-on-a-cluster)
- [Usage Examples](#usage-examples)
- [Troubleshooting](#troubleshooting)

---

## Overview

```
samples.csv
    ↓
┌──────────────────────────────────────────────────┐
│  Input resolution (per sample)                   │
│                                                  │
│  FASTQ ──────────────────────────────────────┐   │
│  BAM ──→ cigarmath/bam2fastx ───────────────-┘   │
└──────────────────────────────────────────────────┘
    ↓
(optional) filtlong downsampling
    ↓
strainline/strainline  →  strainline/{sample_name}/haplotypes.fa
    ↓
strainline/clipqs  →  clipped/{sample_name}.fa
    ↓
(optional) VADR annotation  →  vadr/{sample_name}/
    ↓
(optional) hiv/vadr-tbl2gbk  →  gbk/{sample_name}.gbf + .sqn
    ↓
(optional, per tree_group)
    concat_group_fastas  →  msa/{tree_group}.combined.fa
    MSA/muscle           →  msa/{tree_group}.fa
    phylo/FastTree       →  trees/{tree_group}.nwk
    phylo/phytreeviz     →  figures/{tree_group}.png
```

**Key features:**

- Two flexible input modes: FASTQ or BAM (auto-converted to FASTQ)
- Optional read downsampling with filtlong before assembly
- Haplotype reconstruction via Strainline (long-read-aware, no reference required)
- Reference-based clipping and orientation of haplotypes using clipqs/minimap2
- Three VADR annotation modes: pre-built model, NCBI fetch, or fully offline local files
- Optional cohort-level phylogenetic trees: MSA → FastTree → phytreeviz
  visualization, grouped by the `tree_group` column in `samples.csv`
- Wrappers fetched from GitHub by default; override with a local path for
  development or reproducibility pinning

---

## Pipeline Stages

### 1. Read Preparation

| Rule | Wrapper | When used |
|------|---------|-----------|
| `bam_to_fastq` | `cigarmath/bam2fastx` | `path` column is a `.bam` file |
| *(none)* | — | `path` column is a FASTQ file |

Intermediate FASTQ files from BAM conversion are written to `fastq/`.

### 2. (Optional) Downsampling (`rule downsample`)

Triggered when `downsample_size` is set in config.  Uses the
[snakemake-wrappers `bio/filtlong`](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/filtlong.html)
wrapper with `target_bases = downsample_size`.

If `filtlong_reference` is also set, it is passed as `--assembly` so filtlong
can preferentially retain reads that match the reference genome.

Output: `fastq/{sample_name}.downsampled.fastq`

### 3. Strainline Haplotype Reconstruction (`rule strainline_reconstruct`)

Runs [Strainline](https://github.com/HaploKit/Strainline) on the prepared
reads using the `strainline/strainline` wrapper.

- Requires `strainline_prefix` in config pointing to the Strainline
  installation directory
- 8 threads per sample by default
- Platform: `ont` (Oxford Nanopore)
- Output: `strainline/{sample_name}/haplotypes.fa`

### 4. Strainline Clipping (`rule strainline_clip`)

Aligns and trims haplotypes against the reference FASTA using
`strainline/clipqs`.

- Unaligned terminal ends are clipped
- Minus-strand haplotypes are reverse-complemented
- Haplotypes with insufficient reference coverage (< 20 %) are discarded
- Output: `clipped/{sample_name}.fa`

### 5. (Optional) VADR Annotation

Three sub-steps, each conditional on the annotation mode selected in config
(see [Annotation Modes](#annotation-modes)):

| Rule | Wrapper | Triggered by |
|------|---------|--------------|
| `vadr_genbank` | `hiv/vadr-genbank` | Mode C only |
| `vadr_build` | `hiv/vadr-build` | Modes B and C |
| `vadr_annotate` | `hiv/vadr-annotate` | All modes |

`vadr_annotate` uses `mode='hiv'` which activates `--alt_pass` for all
recoverable HIV-specific alert codes so that proviral sequences with
frameshifts, internal stops, etc. are still classified as passing.

Output: `vadr/{sample_name}/` (per-sample annotation directory)

### 5b. (Optional) GenBank file generation (`rule vadr_tbl2gbk`)

Runs immediately after `vadr_annotate` whenever annotation is enabled.
Uses the `hiv/vadr-tbl2gbk` wrapper to convert the passing sequences FASTA
and NCBI 5-column feature table into an annotated GenBank flat file via
`table2asn`.

- Output: `gbk/{sample_name}.gbf` (GenBank flat file)
- Output: `gbk/{sample_name}.sqn` (ASN.1 submission file for GenBank deposit)
- Set `genbank_organism` in config for proper source feature metadata

### 6. (Optional) MSA per tree_group (`rule msa_by_group`)

Triggered when the `tree_group` column is present in `samples.csv`.
Haplotypes from all samples sharing a `tree_group` label are concatenated
(with sample-name prefixes for traceability) and aligned using
[MUSCLE](https://github.com/rcedgar/muscle) via the `MSA/muscle` wrapper.

- Concatenation: `msa/{tree_group}.combined.fa`
- Aligned MSA: `msa/{tree_group}.fa`

### 7. (Optional) FastTree (`rule fasttree_by_group`)

Infers a maximum-likelihood phylogeny per group using `phylo/FastTree`.

- GTR + Gamma model (appropriate for HIV/retroviral sequences)
- Output: `trees/{tree_group}.nwk` (Newick format)

### 8. (Optional) Tree Visualization (`rule phytreeviz_by_group`)

Visualizes the inferred tree using `phylo/phytreeviz`.

- 300 DPI PNG, 12 × 10 inches
- Branch support values shown
- Output: `figures/{tree_group}.png`

---

## Quick Start

### 1. Prepare your run directory

```bash
mkdir my_reconstruction_run
cd my_reconstruction_run
```

### 2. Create `samples.csv`

**Minimal — FASTQ input, no tree grouping:**

```csv
sample_name,path
patient_01,reads/patient_01.fastq.gz
patient_02,reads/patient_02.fastq.gz
patient_03,reads/patient_03.fastq.gz
```

**BAM input with tree groups:**

```csv
sample_name,path,tree_group
patient_01_wk4,bams/patient_01_wk4.bam,patient_01
patient_01_wk24,bams/patient_01_wk24.bam,patient_01
patient_02_wk4,bams/patient_02_wk4.bam,patient_02
patient_02_wk24,bams/patient_02_wk24.bam,patient_02
```

### 3. Create `run.meta.yaml`

```yaml
strainline_prefix: /data/tools/strainline
reference_fasta:   /data/references/hxb2.fa

# Optional — uncomment to enable
# downsample_size:    500000
# filtlong_reference: /data/references/hxb2.fa

# Annotation — choose one mode (see Annotation Modes section)
# annotation_genbank: /data/vadr_models/K03455      # Mode A: pre-built dir
# annotation_accession: K03455                       # Mode B: fetch from NCBI
# annotation_accession: K03455                       # Mode C: local files
# annotation_local_genbank: refs/K03455.gb
# annotation_local_fasta:   refs/K03455.fasta

samples_csv: samples.csv
# damlab_prefix: /path/to/local/damlab-wrappers
```

### 4. Run the pipeline

```bash
snakemake -s /path/to/damlab-wrappers/workflows/proviral_reconstruction.smk \
          --use-conda --cores 16
```

---

## samples.csv Reference

Every row is one sample.

| Column | Required | Description |
|--------|----------|-------------|
| `sample_name` | yes | Unique name for this sample. Used in all output paths. |
| `path` | yes | Path to a FASTQ (`.fastq` / `.fastq.gz`) or BAM (`.bam`) file. BAM files are automatically converted to FASTQ. |
| `tree_group` | no | Cohort label. Samples sharing the same label are grouped for MSA and phylogenetic tree construction. Omit or leave blank to skip tree-building for that sample. |

---

## Configuration

The pipeline reads `run.meta.yaml` from the working directory by default.
Override with `--configfile` on the command line.

| Key | Required | Default | Description |
|-----|----------|---------|-------------|
| `strainline_prefix` | **yes** | — | Path to the Strainline installation directory. |
| `reference_fasta` | **yes** | — | Reference FASTA used by clipqs for haplotype alignment and orientation. |
| `downsample_size` | no | *(disabled)* | Filtlong `--target_bases` value. When set, reads are downsampled before Strainline. |
| `filtlong_reference` | no | *(disabled)* | FASTA reference for filtlong `--assembly` scoring. Only used when `downsample_size` is also set. |
| `annotation_genbank` | no | *(disabled)* | **Mode A** — path to a pre-built VADR model directory. |
| `annotation_accession` | no | *(disabled)* | **Modes B/C** — bare GenBank accession for model building (e.g. `K03455`). |
| `annotation_local_genbank` | no | *(disabled)* | **Mode C** — path to a local GenBank flat file (`.gb`). Requires `annotation_accession`. |
| `annotation_local_fasta` | no | *(disabled)* | **Mode C** — path to a local FASTA file. Requires `annotation_accession`. |
| `annotation_group` | no | `""` | VADR `--group` label (e.g. `"HIV"`). Used with Modes B and C. |
| `annotation_subgroup` | no | `""` | VADR `--subgroup` label (e.g. `"HIV-1"`). Used with Modes B and C. |
| `genbank_organism` | no | `""` | Organism string for GenBank source records (e.g. `"Human immunodeficiency virus 1"`). |
| `samples_csv` | no | `samples.csv` | Path to the samples CSV relative to the working directory. |
| `damlab_prefix` | no | GitHub `main` | Base URL or local path for damlab-wrappers. |

### `damlab_prefix`

```yaml
# Default — fetch wrappers from GitHub (no entry needed)

# Pin to a specific release tag:
damlab_prefix: https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/tags/v1.2.3

# Use a local checkout during development:
damlab_prefix: /home/user/repos/damlab-wrappers
```

---

## Annotation Modes

### Mode A — Pre-built model directory

Set `annotation_genbank` to the output directory of a previous `vadr-build`
run.  No build step is run; the model directory is used directly.

```yaml
annotation_genbank: /data/vadr_models/K03455
genbank_organism:   Human immunodeficiency virus 1
```

**DAG:** `vadr_annotate` → `vadr_tbl2gbk`

### Mode B — Fetch from NCBI

Set `annotation_accession` to the bare accession (no `.N` version suffix).
`v-build.pl` fetches the sequence and feature table from NCBI at pipeline
start.  **Requires internet access.**

```yaml
annotation_accession: K03455
annotation_group:     HIV
annotation_subgroup:  HIV-1
genbank_organism:     Human immunodeficiency virus 1
```

**DAG:** `vadr_build` → `vadr_annotate` → `vadr_tbl2gbk`

### Mode C — Local files (fully offline)

Provide both a local GenBank flat file and a local FASTA along with
`annotation_accession` (used as the model name).  The `vadr_genbank` rule
normalises the files before building:

- Replaces the LOCUS name with the bare accession (e.g. `HIVHXB2CG` → `K03455`)
- Strips the `.N` version suffix from the VERSION line (`K03455.1` → `K03455`)
- Copies the FASTA unchanged (VADR requires the `>accession.N` header format)

```yaml
annotation_accession:     K03455
annotation_local_genbank: refs/K03455.gb
annotation_local_fasta:   refs/K03455.fasta
annotation_group:         HIV
annotation_subgroup:      HIV-1
genbank_organism:         Human immunodeficiency virus 1
```

**DAG:** `vadr_genbank` → `vadr_build` → `vadr_annotate` → `vadr_tbl2gbk`

> **Custom references with inserted genes:** Mode C is the right choice when
> annotating against a custom reference sequence (e.g. an HIV backbone with a
> reporter gene inserted).  Create a GenBank file for your modified reference
> with the insert annotated as a `CDS` or `misc_feature`, set `LOCUS` and
> `VERSION` to your chosen model name (no `.N` suffix needed for custom files),
> and point `annotation_local_genbank` / `annotation_local_fasta` at those
> files.

---

## Phylogenetic Trees

Tree generation is enabled by adding a `tree_group` column to `samples.csv`.

- All samples **without** a `tree_group` value produce clipped haplotypes
  (and optionally VADR annotations) but are excluded from tree-building.
- A separate tree is built for **each unique `tree_group`** value.
- Haplotype headers in the alignment are prefixed with `{sample_name}__` so
  individual sequences can be identified in the resulting tree.

Example: comparing viral populations at two timepoints across multiple patients:

```csv
sample_name,path,tree_group
p01_wk04,reads/p01_wk04.fastq,cohort_A
p01_wk24,reads/p01_wk24.fastq,cohort_A
p02_wk04,reads/p02_wk04.fastq,cohort_A
p02_wk24,reads/p02_wk24.fastq,cohort_A
```

This produces:

```
msa/cohort_A.fa
trees/cohort_A.nwk
figures/cohort_A.png
```

---

## Output Files

```
fastq/                               # intermediate FASTQs (BAM conversion / downsampling)
  {sample_name}.fastq                #   BAM-to-FASTQ conversion
  {sample_name}.downsampled.fastq    #   filtlong output (if downsample_size set)

strainline/                          # raw Strainline output
  {sample_name}/
    haplotypes.fa                    #   assembled haplotype sequences

clipped/                             # final per-sample haplotypes (always produced)
  {sample_name}.fa                   #   clipped and oriented haplotypes

vadr_model/                          # VADR model files (if annotation is enabled)
  normalised/                        #   Mode C only: normalised GB/FASTA staging area
    {accession}/
      {accession}.gb
      {accession}.fasta
  {accession}/                       #   built model (Modes B and C)
    {accession}.cm
    {accession}.minfo
    ...

vadr/                                # VADR annotation (if annotation enabled)
  {sample_name}/
    *.vadr.pass.fa                   #   sequences passing all VADR checks
    *.vadr.fail.fa                   #   sequences failing one or more checks
    *.vadr.pass.tbl                  #   NCBI 5-column feature table (passing)
    *.vadr.fail.tbl                  #   NCBI 5-column feature table (failing)
    *.vadr.sqc                       #   per-sequence classification summary
    *.vadr.ftr                       #   per-feature annotation table
    *.vadr.alt                       #   alert detail listing

gbk/                                 # GenBank files (if annotation enabled)
  {sample_name}.gbf                  #   annotated GenBank flat file
  {sample_name}.sqn                  #   ASN.1 submission file for NCBI deposit

msa/                                 # multiple sequence alignments (if tree_group set)
  {tree_group}.combined.fa           #   concatenated unaligned haplotypes
  {tree_group}.fa                    #   MUSCLE-aligned MSA

trees/                               # phylogenetic trees (if tree_group set)
  {tree_group}.nwk                   #   Newick tree from FastTree (GTR+Gamma)

figures/                             # tree visualizations (if tree_group set)
  {tree_group}.png                   #   phytreeviz PNG (300 DPI)

logs/                                # all step logs
  {sample_name}.bam2fastq.log
  {sample_name}.filtlong.log
  {sample_name}.strainline.log
  {sample_name}.clipqs.log
  vadr_genbank.log                   #   Mode C normalisation (if applicable)
  vadr_build.log                     #   model build (Modes B and C)
  {sample_name}.vadr_annotate.log
  {sample_name}.tbl2gbk.log
  {tree_group}.muscle.log
  {tree_group}.fasttree.log
  {tree_group}.phytreeviz.log
```

---

## Running on a Cluster

**SLURM example:**

```bash
snakemake -s /path/to/damlab-wrappers/workflows/proviral_reconstruction.smk \
          --use-conda \
          --executor slurm \
          --default-resources slurm_partition=gpu mem_mb=32000 \
          --jobs 20
```

**Resource guidelines:**

| Step | Threads | Memory |
|------|---------|--------|
| `bam_to_fastq` | 1 | 4 GB |
| `downsample` | 1 | 4 GB |
| `strainline_reconstruct` | 8 | 32 GB |
| `strainline_clip` | 1 | 4 GB |
| `vadr_genbank` | 1 | 1 GB |
| `vadr_build` | 1 | 8 GB |
| `vadr_annotate` | 4 | 16 GB |
| `msa_by_group` | 4 | 8 GB |
| `fasttree_by_group` | 4 | 8 GB |
| `phytreeviz_by_group` | 1 | 4 GB |

---

## Usage Examples

### Reconstruction only (no annotation, no trees)

```yaml
# run.meta.yaml
strainline_prefix: /data/tools/strainline
reference_fasta:   /data/references/hxb2.fa
```

```csv
sample_name,path
p01,reads/p01.fastq.gz
p02,reads/p02.fastq.gz
```

### Full pipeline — NCBI model build (Mode B)

```yaml
# run.meta.yaml
strainline_prefix:    /data/tools/strainline
reference_fasta:      /data/references/hxb2.fa
downsample_size:      500000
filtlong_reference:   /data/references/hxb2.fa
annotation_accession: K03455
annotation_group:     HIV
annotation_subgroup:  HIV-1
genbank_organism:     Human immunodeficiency virus 1
```

```csv
sample_name,path,tree_group
p01_early,reads/p01_early.fastq.gz,patient_01
p01_late,reads/p01_late.fastq.gz,patient_01
p02_early,reads/p02_early.fastq.gz,patient_02
p02_late,reads/p02_late.fastq.gz,patient_02
```

### Full pipeline — local files, fully offline (Mode C)

```yaml
# run.meta.yaml
strainline_prefix:        /data/tools/strainline
reference_fasta:          /data/references/hxb2.fa
annotation_accession:     K03455
annotation_local_genbank: /data/refs/K03455.gb
annotation_local_fasta:   /data/refs/K03455.fasta
annotation_group:         HIV
annotation_subgroup:      HIV-1
genbank_organism:         Human immunodeficiency virus 1
```

### Custom reference with an inserted gene (Mode C)

For a modified HIV sequence with a reporter gene (or any custom reference),
create a GenBank file for your reference with the insert annotated.  Use a
simple model name with no `.N` suffix:

```
LOCUS       MY_HIV_GFP   9900 bp
...
VERSION     MY_HIV_GFP
```

```yaml
annotation_accession:     MY_HIV_GFP
annotation_local_genbank: refs/my_hiv_gfp.gb
annotation_local_fasta:   refs/my_hiv_gfp.fa
```

### Re-use a previously built model (Mode A)

```yaml
# After running Mode B or C once, skip future builds with:
annotation_genbank: vadr_model/K03455
```

---

## Troubleshooting

### `strainline_prefix` not found

```
ValueError: Config key 'strainline_prefix' is required.
```

Set `strainline_prefix` in `run.meta.yaml` to the directory where Strainline
is installed, e.g. `/data/tools/strainline`.  Strainline is not conda-installable
and must be installed manually.

### Strainline produces no haplotypes

This can happen with very low read depth.  Try reducing `downsample_size` or
removing it entirely to give Strainline the maximum read count.

### clipqs removes all sequences

The `clipped/{sample_name}.fa` output is empty or missing.  This usually
means the reference FASTA (`reference_fasta`) does not match your data
(e.g. wrong subtype).  Verify orientation by mapping your reads to the
reference with minimap2 and inspecting the alignment.

### VADR annotation directory already exists

```
FileExistsError: Output directory 'vadr/sample_name' already exists and is non-empty.
```

`v-annotate.pl` must create its output directory from scratch.  Remove the
existing directory and re-run:

```bash
rm -rf vadr/{sample_name}
snakemake --use-conda --cores 4 vadr/{sample_name}
```

### Mode C — VADR version/accession mismatch

If `vadr_build` still fails after `vadr_genbank` with
`version/accession mismatch`, check that `annotation_accession` matches the
`ACCESSION` line (not the `LOCUS` name) in your GenBank file.  The
`vadr_genbank` wrapper will replace the LOCUS name and strip the VERSION
suffix but it cannot fix a file where the accession field itself is wrong.

### Wrappers not found (local `damlab_prefix`)

Ensure `damlab_prefix` points to the **root** of the `damlab-wrappers`
checkout, not to a sub-directory:

```yaml
damlab_prefix: /home/user/repos/damlab-wrappers   # correct
damlab_prefix: /home/user/repos/damlab-wrappers/workflows  # wrong
```
