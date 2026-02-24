# Proviral NFL Pipeline

A comprehensive Snakemake pipeline for processing Nanopore sequencing data from Oxford Nanopore Technologies (ONT) sequencers.
This pipeline handles basecalling, demultiplexing, alignment, and quality control for viral near-full-length (NFL) sequencing projects.

## Table of Contents

- [Overview](#overview)
- [Pipeline Stages](#pipeline-stages)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Profile Configuration](#profile-configuration)
- [Input Modes](#input-modes)
- [Basecalling Modes](#basecalling-modes)
- [Output Files](#output-files)
- [Usage Examples](#usage-examples)
- [Troubleshooting](#troubleshooting)

## Overview

The proviral NFL pipeline processes Nanopore sequencing data through the following stages:

```
POD5 Files → Basecalling → Demultiplexing → Read Groups → Alignment → Sorting → Analysis → Metrics → Report
```

**Key Features:**
- Supports duplex and simplex basecalling with Dorado
- Scatter mode for distributed cluster execution
- Automatic demultiplexing with multiple barcode kits
- Reference-based alignment with minimap2
- Haplotype reconstruction with Strainline
- Deletion block detection for identifying large deletions
- Comprehensive QC metrics with MultiQC
- Flexible configuration for different compute environments

## Pipeline Stages

### 1. Basecalling (`rules/basecalling.smk`)

Converts raw POD5 signal data to DNA sequences using Oxford Nanopore's Dorado basecaller.

**Modes:**
- **Duplex mode** (default): Higher accuracy using duplex reads
- **Simplex mode**: Faster basecalling with standard accuracy

**Execution Modes:**
- **Non-scattered** (default): All POD5 files processed in single GPU job
- **Scattered**: Each POD5 file processed independently (optimal for SLURM clusters)

**Rules:**
- `pod5_duplexing` / `pod5_simplex` - Standard basecalling
- `pod5_duplexing_scattered` / `pod5_simplex_scattered` - Scattered basecalling
- `merge_duplex_scattered` / `merge_simplex_scattered` - Merge scattered outputs

### 2. Demultiplexing (`rules/demux.smk`)

Separates multiplexed samples based on barcode sequences.

**Supported Kits:**
- SQK-NBD114-24 (Native Barcoding Kit 24)
- SQK-NBD111-24 (Native Barcoding Kit 24 v11)
- EXP-NBD103, EXP-NBD104, EXP-NBD114
- SQK-RBK001, SQK-RBK004

**Rule:**
- `demux_run` - Demultiplexes basecalled reads into per-sample BAM files

### 3. Alignment (`rules/alignment.smk`)

Aligns demultiplexed reads to reference sequences and adds read group information.

**Rules:**
- `add_read_groups` - Adds sample metadata to BAM files
- `align_reads` - Aligns reads using minimap2 via Dorado aligner
- `sort_reads` - Sorts aligned BAM files
- `index_reads` - Creates BAM index files

### 4. Metrics (`rules/metrics.smk`)

Generates quality control metrics for each sample.

**Rules:**
- `metrics` - Calculates HIV-specific metrics (can be adapted for other viruses)
- `depth` - Computes per-position coverage depth
- `samtools_stats` - Standard BAM statistics

### 5. Analysis (`rules/analysis.smk`)

Performs advanced analysis on aligned reads.

**Rules:**
- `bam_to_fastq` - Converts aligned BAM to FASTQ for strainline input
- `strainline` - Haplotype reconstruction using Strainline (identifies viral quasispecies)
- `deletion_block_detection` - Identifies large deletion blocks in aligned reads (useful for detecting defective proviruses)

### 6. Reporting (`rules/reporting.smk`)

Aggregates all QC metrics into interactive HTML report.

**Rules:**
- `run_report` - Generates MultiQC report combining all metrics

## Quick Start

### Minimal Setup

1. **Prepare your directory:**
   ```bash
   mkdir my_run
   cd my_run
   ```

2. **Create run configuration** (`run.meta.yaml`):
   ```yaml
   run_date: "2026-02-03"
   flow_cell_id: "FAT12345"
   flow_cell_product_code: "FLO-MIN114"
   position: "MN35859"
   run_id: "my_run_001"
   project_code: "NEW"
   
   samples_csv: samples.csv
   library_prep_kits:
     - SQK-NBD114-24
   
   DORADO_MODE: duplex
   DORADO_MODEL: sup
   MODEL_ROOT: /path/to/dorado_models/
   REFERENCE: /path/to/reference.fasta
   pod5_path: /path/to/pod5/
   ```

3. **Create samples file** (`samples.csv`):
   ```csv
   sample_name,barcode,reference
   sample1,barcode01,ref.fasta
   sample2,barcode02,ref.fasta
   sample3,barcode03,ref.fasta
   ```

4. **Run pipeline:**
   ```bash
   snakemake --workflow proviral_nfl.smk --cores 8 --use-conda
   ```

## Configuration

### Required Fields

The pipeline uses a YAML configuration file (default: `run.meta.yaml`) with the following required fields:

```yaml
# Run metadata
run_date: "YYYY-MM-DD"           # Date of sequencing run
flow_cell_id: "FAT12345"         # Flowcell identifier
flow_cell_product_code: "FLO-MIN114"  # One of: FLO-MIN106, FLO-MIN111, FLO-MIN114
position: "MN35859"              # MinION position ID
run_id: "unique_run_id"          # Unique identifier for this run
project_code: "PROJECT"          # Project code for billing/tracking

# Sample information
samples_csv: samples.csv         # Path to samples CSV file
library_prep_kits:               # Barcode kits used
  - SQK-NBD114-24

# Basecalling parameters
DORADO_MODE: duplex              # duplex or simplex
DORADO_MODEL: sup                # Model accuracy: fast, hac, sup
MODEL_ROOT: /path/to/models/     # Dorado models directory
REFERENCE: /path/to/ref.fasta    # Reference sequence for alignment
pod5_path: /path/to/pod5/        # Directory containing POD5 files
```

### Optional Fields

```yaml
# Operator/technician information (optional)
flow_cell_operator_initials: "ABC"
library_prep_initials: "XYZ"
pcr_initials: "DEF"
data_processing_initials: "GHI"

# PCR protocol (optional)
pcr_protocol: "hiv-nfl"          # One of: hiv-nfl, hiv-nfl-umi, hiv-targeted, siv-nfl, siv-nfl-umi

# Advanced options
snakemake_wrapper_tag: v8.1.1    # Snakemake-wrappers version
damlab_prefix: /path/to/wrappers # Path to damlab-wrappers directory

# Scatter mode (for clusters)
DORADO_SCATTER: False            # Enable scatter basecalling
SCATTER_THREADS: 12              # Threads per scattered job
MERGE_THREADS: 16                # Threads for merge step

# Alternative input mode
demuxed_bam_path: /path/to/bam   # Skip basecalling, use pre-demuxed BAM

# Analysis options
STRAINLINE_PREFIX: /path/to/strainline  # Path to strainline installation
MIN_DELETION_SIZE: 50                    # Minimum deletion size to detect (default: 50)
```

### Samples CSV Format

The `samples.csv` file defines sample-to-barcode mapping:

```csv
sample_name,barcode,reference
patient_001_baseline,barcode01,HIV_HXB2.fasta
patient_001_week4,barcode02,HIV_HXB2.fasta
patient_002_baseline,barcode03,HIV_HXB2.fasta
```

**Columns:**
- `sample_name` - Unique sample identifier (used in output filenames)
- `barcode` - Barcode ID matching the kit (e.g., barcode01-barcode24)
- `reference` - Reference genome for alignment (can be same for all samples)

## Profile Configuration

Profiles define compute resource settings for different execution environments.

### Local Execution Profile

For running on a single workstation:

```yaml
# profiles/local/config.yaml
jobs: 4
use-conda: True

config:
  MODEL_ROOT: /opt/dorado_models/
  DORADO_SCATTER: False  # Single GPU, no scatter needed

default-resources:
  mem_mb: 8000
  runtime: 240  # 4 hours
```

### SLURM Cluster Profile (Non-Scatter)

For cluster execution without scatter mode:

```yaml
# profiles/cluster/config.yaml
executor: slurm
jobs: 32
use-conda: True
conda-prefix: /shared/conda_envs/

config:
  MODEL_ROOT: /shared/dorado_models/
  DORADO_SCATTER: False

default-resources:
  slurm_account: myproject
  runtime: 240

set-resources:
  pod5_duplexing:
    mem_mb: 16000
    runtime: 1440  # 24 hours
    nodes: 1
    cpus_per_task: 4
    slurm_partition: 'gpu'
    slurm_extra: "--gres=gpu:2"
```

### SLURM Cluster Profile (Scatter Mode)

For distributed basecalling across multiple nodes:

```yaml
# profiles/cluster_scatter/config.yaml
executor: slurm
jobs: 100
use-conda: True
conda-prefix: /shared/conda_envs/

config:
  MODEL_ROOT: /shared/dorado_models/
  DORADO_SCATTER: True
  SCATTER_THREADS: 1
  MERGE_THREADS: 16

default-resources:
  slurm_account: myproject
  runtime: 240

set-resources:
  # Scattered basecalling rules - each POD5 file gets own GPU job
  pod5_duplexing_scattered:
    mem_mb: 4000
    runtime: 1440  # 24 hours
    nodes: 1
    cpus_per_task: 1
    slurm_partition: 'gpu'
    slurm_extra: "--gres=gpu:1"
  
  pod5_simplex_scattered:
    mem_mb: 4000
    runtime: 1440
    nodes: 1
    cpus_per_task: 1
    slurm_partition: 'gpu'
    slurm_extra: "--gres=gpu:1"
  
  # Merge rules - CPU only, no GPU needed
  merge_duplex_scattered:
    mem_mb: 8000
    runtime: 120  # 2 hours
    nodes: 1
    cpus_per_task: 16
    slurm_partition: 'standard'
  
  merge_simplex_scattered:
    mem_mb: 8000
    runtime: 120
    nodes: 1
    cpus_per_task: 16
    slurm_partition: 'standard'
  
  # Other rules
  demux_run:
    mem_mb: 16000
    runtime: 240
    cpus_per_task: 16
  
  align_reads:
    mem_mb: 8000
    runtime: 120
    cpus_per_task: 4
```

## Input Modes

The pipeline supports two input modes, automatically detected:

### 1. POD5 Mode (Default)

Start from raw POD5 signal files. The pipeline will basecall and demultiplex.

**Configuration:**
```yaml
pod5_path: /path/to/pod5_directory/
DORADO_MODE: duplex
DORADO_MODEL: sup
```

**Directory structure:**
```
pod5_directory/
├── channel-1.pod5
├── channel-2.pod5
└── channel-N.pod5
```

### 2. Pre-Demuxed BAM Mode

Skip basecalling/demultiplexing if you already have demuxed BAM files.

**Configuration:**
```yaml
demuxed_bam_path: /path/to/demuxed_bams/
```

**Directory structure:**
```
demuxed_bams/
├── sample1.bam
├── sample2.bam
└── sample3.bam
```

## Basecalling Modes

### Duplex vs Simplex

**Duplex Mode** (recommended for accuracy):
```yaml
DORADO_MODE: duplex
DORADO_MODEL: sup  # or: fast, hac
```
- Higher accuracy using complementary strand information
- Slower processing time
- Best for final analysis

**Simplex Mode** (faster):
```yaml
DORADO_MODE: simplex
DORADO_MODEL: sup  # or: fast, hac
```
- Standard basecalling accuracy
- Faster processing
- Good for QC or preliminary analysis

### Scatter Mode

Scatter mode distributes basecalling across multiple GPU nodes for faster processing.

#### When to Use Scatter Mode

**Use Scatter Mode When:**
- Running on SLURM or similar cluster
- Have pre-split POD5 files (by channel)
- Need faster turnaround for large datasets
- Have many GPU nodes available

**Use Non-Scatter Mode When:**
- Running on single GPU workstation
- Have few POD5 files
- Prefer simplicity over speed

#### Enabling Scatter Mode

**Step 1: Pre-split POD5 files**

Split your POD5 files by channel for optimal parallelization:

```bash
# Using pod5 command line tools
pod5 view -r input.pod5 --include 'read_id,channel' --output read2channel.tsv
pod5 subset -r input.pod5 --summary read2channel.tsv --columns channel --output split_pod5/
```

Or using the pipeline wrapper (create a separate Snakefile):

```python
rule split_pod5:
    input:
        pod5_dir="raw_pod5/"
    output:
        directory("split_pod5/")
    wrapper:
        "file://path/to/damlab-wrappers/pod5/split_by_channel"
```

**Step 2: Enable in configuration**

```yaml
pod5_path: /path/to/split_pod5/  # Directory with split files
DORADO_SCATTER: True
SCATTER_THREADS: 1                # Threads per scattered job
MERGE_THREADS: 16                 # Threads for final merge
```

**Step 3: Use appropriate profile**

```bash
snakemake --profile profiles/cluster_scatter/
```

#### Scatter Mode Pipeline Flow

```
Non-Scatter:
  pod5_duplexing (all files, 1 job) → basecalled.bam → demux → ...

Scatter:
  pod5_discovery (checkpoint) →
    pod5_duplexing_scattered (file 1, job 1) → scattered/channel-1.bam
    pod5_duplexing_scattered (file 2, job 2) → scattered/channel-2.bam
    pod5_duplexing_scattered (file N, job N) → scattered/channel-N.bam
  → merge_duplex_scattered (all scattered BAMs) → basecalled.bam → demux → ...
```

**Performance Benefits:**
- Linear speedup with number of GPU nodes
- Each POD5 file scheduled independently
- Better cluster resource utilization
- Can prioritize/cancel individual file jobs

## Output Files

The pipeline generates the following directory structure:

```
output/
├── duplex/                       # or simplex/
│   ├── basecalled.bam            # Basecalled reads (all barcodes)
│   └── scattered/                # Intermediate scattered BAMs (scatter mode only)
│       ├── channel-1.bam
│       └── channel-N.bam
├── demux/                        # Demultiplexed reads per sample
│   ├── sample1.bam
│   ├── sample2.bam
│   └── sample3.bam
├── labeled/                      # BAMs with read groups added
│   ├── sample1.rg.bam
│   └── ...
├── aligned/                      # Aligned and sorted BAMs
│   ├── sample1.sorted.bam
│   ├── sample1.sorted.bam.bai
│   ├── sample1.mm2.bam           # Intermediate unsorted BAM
│   └── ...
├── metrics/                      # Per-sample QC metrics
│   ├── sample1.hivmetrics.yaml
│   ├── sample1.depth.txt
│   └── ...
├── samtools_stats/               # SAMtools statistics
│   ├── sample1.txt
│   └── ...
├── analysis/                     # Advanced analysis outputs
│   ├── sample1.haplotypes.fa     # Reconstructed haplotypes (Strainline)
│   ├── sample1.deletion_reads.csv    # Per-read deletion info
│   ├── sample1.deletion_blocks.csv   # Unique deletion blocks found
│   ├── sample1.deletion_summary.yaml # Summary statistics
│   └── ...
└── qc/                           # Final QC report
    ├── multiqc.html              # Interactive QC report
    └── multiqc_data/             # Supporting data files
```

### Key Output Files

**Primary Outputs:**
- `aligned/{sample}.sorted.bam` - Final aligned BAM files for analysis
- `analysis/{sample}.haplotypes.fa` - Reconstructed viral haplotypes
- `analysis/{sample}.deletion_blocks.csv` - Detected deletion blocks
- `qc/multiqc.html` - Comprehensive QC report (open in browser)

**Intermediate Files:**
- `duplex/basecalled.bam` - All basecalled reads before demux
- `demux/{sample}.bam` - Demultiplexed reads per sample
- `metrics/{sample}.depth.txt` - Per-position coverage
- `analysis/{sample}.deletion_summary.yaml` - Deletion statistics for MultiQC

## Usage Examples

### Example 1: Local Workstation (Single GPU)

```bash
# Directory setup
cd /data/my_experiment
mkdir -p pod5/

# Configuration
cat > run.meta.yaml << EOF
run_date: "2026-02-03"
flow_cell_id: "FAT12345"
flow_cell_product_code: "FLO-MIN114"
position: "MN35859"
run_id: "exp001"
project_code: "NEW"

samples_csv: samples.csv
library_prep_kits:
  - SQK-NBD114-24

DORADO_MODE: duplex
DORADO_MODEL: sup
MODEL_ROOT: /opt/dorado_models/
REFERENCE: /data/refs/HIV_HXB2.fasta
pod5_path: pod5/
EOF

# Run
snakemake --snakefile /path/to/damlab-wrappers/workflows/proviral_nfl.smk \
  --cores 8 --use-conda
```

### Example 2: SLURM Cluster (Scatter Mode)

```bash
# Pre-split POD5 files
pod5 view -r raw.pod5 --include 'read_id,channel' --output read2channel.tsv
pod5 subset -r raw.pod5 --summary read2channel.tsv --columns channel --output split_pod5/

# Create run config (run.meta.yaml)
cat > run.meta.yaml << EOF
run_date: "2026-02-03"
flow_cell_id: "FAT67890"
flow_cell_product_code: "FLO-MIN114"
position: "MN35859"
run_id: "exp002"
project_code: "MyProject"

samples_csv: samples.csv
library_prep_kits:
  - SQK-NBD114-24

DORADO_MODE: duplex
DORADO_MODEL: sup
REFERENCE: /shared/refs/HIV_HXB2.fasta
pod5_path: split_pod5/
EOF

# Run with cluster profile (DORADO_SCATTER set in profile)
snakemake --snakefile /path/to/damlab-wrappers/workflows/proviral_nfl.smk \
  --profile /path/to/profiles/cluster_scatter/
```

### Example 3: Pre-Demuxed Input

Skip basecalling if you have pre-processed BAMs:

```bash
cat > run.meta.yaml << EOF
run_date: "2026-02-03"
flow_cell_id: "FAT99999"
flow_cell_product_code: "FLO-MIN114"
position: "MN35859"
run_id: "exp003"
project_code: "NEW"

samples_csv: samples.csv
library_prep_kits:
  - SQK-NBD114-24

REFERENCE: /data/refs/HIV_HXB2.fasta
demuxed_bam_path: /data/preprocessed_bams/  # Skip basecalling
EOF

snakemake --snakefile /path/to/damlab-wrappers/workflows/proviral_nfl.smk \
  --cores 16 --use-conda
```

### Example 4: Dry Run (Check Pipeline)

Always test with dry run first:

```bash
# See what rules will run
snakemake --snakefile /path/to/damlab-wrappers/workflows/proviral_nfl.smk \
  --dry-run --printshellcmds

# Generate rule graph
snakemake --snakefile /path/to/damlab-wrappers/workflows/proviral_nfl.smk \
  --rulegraph | dot -Tpng > pipeline.png

# Generate DAG for specific sample
snakemake --snakefile /path/to/damlab-wrappers/workflows/proviral_nfl.smk \
  --dag aligned/sample1.sorted.bam | dot -Tpng > dag.png
```

## Troubleshooting

### Common Issues

**Issue: "No pod5 files found"**
```
Solution: Check that pod5_path points to correct directory
- Verify: ls /path/to/pod5/*.pod5
- Config should have: pod5_path: /full/path/to/pod5/
```

**Issue: "REFERENCE not specified"**
```
Solution: Add REFERENCE to config
- Example: REFERENCE: /path/to/ref.fasta
- Must be a valid FASTA file
```

**Issue: Scatter mode not working**
```
Solution: Check requirements
1. POD5 files should be pre-split
2. DORADO_SCATTER: True in config or profile
3. Profile has set-resources for scattered rules
4. Check: snakemake --dry-run
```

**Issue: GPU not detected**
```
Solution: Check Dorado installation
- Test: dorado basecaller --device cuda:all --help
- Verify CUDA drivers installed
- Check GPU availability: nvidia-smi
```

**Issue: Barcode kit not recognized**
```
Solution: Use correct kit name
- Valid kits: SQK-NBD114-24, SQK-NBD111-24, EXP-NBD103, etc.
- Check dorado demux --help for full list
```

**Issue: Memory errors during merge**
```
Solution: Increase MERGE_THREADS resources
- In profile: merge_duplex_scattered: mem_mb: 16000
- Or reduce number of input files
```

### Debugging Tips

**Check rule execution:**
```bash
# Verbose output
snakemake --verbose

# Show shell commands
snakemake --printshellcmds

# Keep going despite failures
snakemake --keep-going
```

**Check specific rule:**
```bash
# Force re-run of specific rule
snakemake --forcerun pod5_duplexing

# Target specific output
snakemake aligned/sample1.sorted.bam
```

**Check logs:**
```
Logs are in same directory as outputs:
- demux.log - Demultiplexing log
- aligned/{sample}.log - Alignment logs
- metrics/{sample}.hivmetrics.log - Metrics logs
- analysis/{sample}.strainline.log - Strainline logs
- analysis/{sample}.deletion_detection.log - Deletion detection logs
```

**Validate configuration:**
```bash
# Check config loading
snakemake --configfile run.meta.yaml --printshellcmds --dry-run

# List all rules
snakemake --list

# List all targets
snakemake --list-target-rules
```

## Advanced Topics

### Custom Model Path

If using non-standard Dorado models:

```yaml
MODEL_ROOT: /custom/path/to/models/
DORADO_MODEL: dna_r10.4.1_e8.2_5khz_stereo@v1.3  # Full model name
```

### Multiple References

If different samples need different references, specify per-sample in CSV:

```csv
sample_name,barcode,reference
hiv_patient1,barcode01,/refs/HIV_HXB2.fasta
siv_sample1,barcode02,/refs/SIV_MAC239.fasta
```

### Resource Tuning

Adjust resources based on your data:

```yaml
# In profile config.yaml
set-resources:
  pod5_duplexing:
    runtime: 2880  # 48 hours for very large runs
    mem_mb: 32000  # More memory for complex data
```

### Conda Environment Control

```bash
# Use specific conda environment
snakemake --use-conda --conda-prefix /path/to/envs/

# Create environments before running
snakemake --conda-create-envs-only
```

## Related Documentation

- [Basecalling Configuration Examples](rules/basecalling_example_configs.md)
- [Scatter Mode Implementation](SCATTER_MODE_IMPLEMENTATION.md)
- [Dorado Documentation](https://github.com/nanoporetech/dorado)
- [Snakemake Profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)

## Citation

If you use this pipeline in your research, please cite:
- The damlab-wrappers repository
- Dorado basecaller (Oxford Nanopore Technologies)
- Snakemake workflow management system

## Support

For issues or questions:
- Check the troubleshooting section above
- Review example configurations in `workflows/test_configs/`
- Consult wrapper-specific README files in each wrapper directory
- Validate your setup with `--dry-run` before full execution

---

**Version:** 2.1.0 (with strainline and deletion detection)  
**Last Updated:** February 24, 2026  
**Maintainer:** Dampier Lab
