"""proviral_crispr.smk — CRISPResso automation pipeline.

Automates CRISPResso (and optional CRISPRessoCompare) analysis from flexible
input types: paired or single-end FASTQ, or a BAM file (with optional region
slicing via cigarmath/slice).

Config keys (config.yaml or run.meta.yaml):
    samples_csv     : path to samples CSV (default: samples.csv)
    MIN_DELETION_SIZE       : optional; min deletion length for deletion_block_detection (default: 50)
    DELETION_MERGE_DISTANCE : optional; merge nearby deletion blocks (default: 10)
    damlab_prefix   : base location for damlab-wrappers. Can be:
                        - a local filesystem path  (e.g. /path/to/damlab-wrappers)
                        - a URL                    (e.g. https://raw.githubusercontent.com/...)
                      Default: https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main
                      Local paths are prefixed with 'file:' in the Snakemake wrapper directive;
                      URLs are used as-is.

samples.csv columns:
    sample_name  (required) unique run name
    grna         (required) guide RNA sequence
    amplicon     (required) amplicon sequence string OR path to a FASTA file
    fastq_r1     (cond.)   R1 FASTQ (fastq mode)
    fastq_r2     (optional) R2 FASTQ for paired-end
    bam_file     (cond.)   BAM file (bam mode)
    region       (optional) chr:start-stop, only used with bam_file
    comparison   (optional) 'experiment' or 'control' — enables CRISPRessoCompare
                            every experiment sample is compared to every control
    deletion_query (optional) regions for per-region deletion/coverage stats on the BAM
                            (ref:start-end; multiple separated by ';'). Passed to
                            cigarmath/deletion_block_detection params.query. BAM samples only.
"""

import os
import itertools
from os.path import join

import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

configfile: "run.meta.yaml"

WORKFLOW_DIR = workflow.basedir
_GITHUB_DEFAULT = (
    "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main"
)
DL_PREFIX = config.get("damlab_prefix", _GITHUB_DEFAULT)


def wrapper_path(subpath):
    """Return a Snakemake wrapper path for a given subpath within damlab-wrappers.

    Local paths are prefixed with 'file:'; URLs are used verbatim so Snakemake
    fetches them directly from the remote host.
    """
    if DL_PREFIX.startswith("http://") or DL_PREFIX.startswith("https://"):
        return f"{DL_PREFIX}/{subpath}"
    return f"file:{DL_PREFIX}/{subpath}"

# ---------------------------------------------------------------------------
# Load samples
# ---------------------------------------------------------------------------

SAMPLES = pd.read_csv(config.get("samples_csv", "samples.csv"))

# ---------------------------------------------------------------------------
# Comparison detection (optional column)
# ---------------------------------------------------------------------------

HAS_COMPARISON = (
    "comparison" in SAMPLES.columns
    and SAMPLES["comparison"].notna().any()
)

if HAS_COMPARISON:
    EXPERIMENT_SAMPLES = (
        SAMPLES[SAMPLES["comparison"] == "experiment"]["sample_name"].tolist()
    )
    CONTROL_SAMPLES = (
        SAMPLES[SAMPLES["comparison"] == "control"]["sample_name"].tolist()
    )
else:
    EXPERIMENT_SAMPLES = []
    CONTROL_SAMPLES = []

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_sample(wildcards):
    """Return the row for the given sample_name wildcard."""
    rows = SAMPLES[SAMPLES["sample_name"] == wildcards.sample_name]
    if rows.empty:
        raise ValueError(f"Sample '{wildcards.sample_name}' not found in samples CSV.")
    return rows.iloc[0]


def _notna(value):
    """Return True if value is a non-empty, non-NaN string."""
    try:
        return pd.notna(value) and str(value).strip() != ""
    except Exception:
        return False


if "bam_file" in SAMPLES.columns:
    BAM_SAMPLE_NAMES = SAMPLES.loc[
        SAMPLES["bam_file"].apply(_notna), "sample_name"
    ].tolist()
else:
    BAM_SAMPLE_NAMES = []


def get_input_mode(row):
    """Return 'bam_region', 'bam', or 'fastq' for a sample row."""
    if _notna(row.get("bam_file")):
        if _notna(row.get("region")):
            return "bam_region"
        return "bam"
    return "fastq"


def is_amplicon_file(value):
    """Return True if the amplicon value is a path to an existing file."""
    return os.path.exists(str(value))


def _opt(wildcards, col):
    """Return the value of an optional samples.csv column, or None if absent/NaN."""
    row = get_sample(wildcards)
    v = row.get(col) if hasattr(row, "get") else None
    if not _notna(v):
        return None
    return v


def _opt_deletion_query(wildcards):
    """``deletion_query`` cell as a trimmed string (always ``str()`` for pandas/numpy scalars)."""
    row = get_sample(wildcards)
    v = row.get("deletion_query") if hasattr(row, "get") else None
    if not _notna(v):
        return None
    return str(v).strip()


# ---------------------------------------------------------------------------
# Input functions for the crispresso rule
# ---------------------------------------------------------------------------

def get_r1_fastq(wildcards):
    """Return the R1 FASTQ path for a sample (bam-derived or direct)."""
    row = get_sample(wildcards)
    mode = get_input_mode(row)
    if mode == "bam_region":
        return f"fastq/{wildcards.sample_name}.slice.fastq"
    if mode == "bam":
        return f"fastq/{wildcards.sample_name}.bam.fastq"
    return row["fastq_r1"]


def get_r2_fastq(wildcards):
    """Return R2 FASTQ path, or [] if single-end or BAM-derived."""
    row = get_sample(wildcards)
    mode = get_input_mode(row)
    if mode != "fastq":
        return []
    r2 = row.get("fastq_r2")
    if _notna(r2):
        return [r2]
    return []


def get_amplicon_fasta(wildcards):
    """Return the amplicon FASTA path, or [] if amplicon is a sequence string."""
    row = get_sample(wildcards)
    amplicon = row["amplicon"]
    if is_amplicon_file(amplicon):
        return [str(amplicon)]
    return []


# ---------------------------------------------------------------------------
# Output collection helpers
# ---------------------------------------------------------------------------

def get_all_crispresso_outputs(wildcards):
    return [
        f"crispresso/CRISPResso_on_{s}"
        for s in SAMPLES["sample_name"].tolist()
    ]


def get_all_compare_outputs(wildcards):
    if not HAS_COMPARISON:
        return []
    return [
        f"crispresso/CRISPRessoCompare_{exp}_vs_{ctrl}"
        for exp, ctrl in itertools.product(EXPERIMENT_SAMPLES, CONTROL_SAMPLES)
    ]


def get_all_deletion_outputs(wildcards):
    """Deletion block detection outputs (BAM-backed samples only)."""
    paths = []
    for s in BAM_SAMPLE_NAMES:
        paths.extend(
            [
                f"deletion_detection/{s}.deletion_reads.csv",
                f"deletion_detection/{s}.deletion_blocks.csv",
                f"deletion_detection/{s}.deletion_summary.yaml",
                f"deletion_detection/{s}.deletion_query_stats.csv",
            ]
        )
    return paths


def get_all_outputs(wildcards):
    return (
        get_all_crispresso_outputs(wildcards)
        + get_all_compare_outputs(wildcards)
        + ["crispresso/CRISPRessoAggregate_on_all"]
        + get_all_deletion_outputs(wildcards)
    )


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

rule all:
    input:
        get_all_outputs


rule slice_bam_region:
    """Extract reads overlapping a region from a BAM file (cigarmath/slice).

    Reads are sliced so only bases covering the target region are returned,
    which is appropriate when long reads span a larger context than the amplicon.
    Used when bam_file AND region are both provided in samples.csv.
    """
    input:
        lambda wc: get_sample(wc)["bam_file"],
    output:
        fastq = "fastq/{sample_name}.slice.fastq",
    params:
        region = lambda wc: get_sample(wc)["region"],
        sample_name = lambda wc: wc.sample_name,
    log:
        "logs/{sample_name}.slice.log",
    wrapper:
        wrapper_path("cigarmath/slice")


rule bam_to_fastq:
    """Convert a BAM file to FASTQ (cigarmath/bam2fastx).

    Full read sequences are exported (no region slicing).
    Used when bam_file is provided but no region is specified in samples.csv.
    """
    input:
        lambda wc: get_sample(wc)["bam_file"],
    output:
        fastq = "fastq/{sample_name}.bam.fastq",
    params:
        mapped_only = True,
    log:
        "logs/{sample_name}.bam2fastx.log",
    wrapper:
        wrapper_path("cigarmath/bam2fastx")


rule deletion_block_detection:
    """Detect reference deletion blocks per sample (cigarmath/deletion_block_detection).

    Runs only for rows with ``bam_file`` set. Uses the same BAM as slice/bam2fastq.
    Optional ``deletion_query`` column is passed as ``params.deletion_query`` (regions
    ``ref:start-end``, multiple separated by semicolons); see wrapper README.
    """
    input:
        bams=lambda wc: get_sample(wc)["bam_file"],
    output:
        reads="deletion_detection/{sample_name}.deletion_reads.csv",
        deletions="deletion_detection/{sample_name}.deletion_blocks.csv",
        summary="deletion_detection/{sample_name}.deletion_summary.yaml",
        query_stats="deletion_detection/{sample_name}.deletion_query_stats.csv",
    params:
        min_deletion_size=config.get("MIN_DELETION_SIZE", 50),
        merge_distance=config.get("DELETION_MERGE_DISTANCE", 10),
        sample_name=lambda wc: wc.sample_name,
        deletion_query=lambda wc: _opt_deletion_query(wc),
    log:
        "logs/{sample_name}.deletion_detection.log",
    wrapper:
        wrapper_path("cigarmath/deletion_block_detection")


rule crispresso:
    """Run CRISPResso on a single sample.

    Input reads come from:
      - FASTQ mode:       fastq_r1 (and optional fastq_r2) directly from CSV
      - BAM mode:         fastq/{sample_name}.bam.fastq  (via bam_to_fastq)
      - BAM+region mode:  fastq/{sample_name}.slice.fastq (via slice_bam_region)

    The amplicon is resolved from the 'amplicon' column:
      - If the value is an existing file path, it is passed as input.amplicon_fasta.
      - Otherwise the value is treated as the amplicon sequence string.

    Optional per-sample columns forwarded to the wrapper when present and
    non-empty: guide_name, amplicon_name, quantification_window_center,
    quantification_window_size, expected_hdr_amplicon_seq.
    """
    input:
        fastq_r1       = get_r1_fastq,
        fastq_r2       = get_r2_fastq,
        amplicon_fasta = get_amplicon_fasta,
    output:
        directory("crispresso/CRISPResso_on_{sample_name}"),
    params:
        amplicon_seq                 = lambda wc: None if is_amplicon_file(get_sample(wc)["amplicon"]) else str(get_sample(wc)["amplicon"]),
        guide_seq                    = lambda wc: get_sample(wc)["grna"],
        guide_name                   = lambda wc: _opt(wc, "guide_name"),
        amplicon_name                = lambda wc: _opt(wc, "amplicon_name"),
        quantification_window_center = lambda wc: _opt(wc, "quantification_window_center"),
        quantification_window_size   = lambda wc: _opt(wc, "quantification_window_size"),
        expected_hdr_amplicon_seq    = lambda wc: _opt(wc, "expected_hdr_amplicon_seq"),
    threads: 4
    log:
        "logs/{sample_name}.crispresso.log",
    wrapper:
        wrapper_path("CRISPR/crispresso-core")


rule crispresso_compare:
    """Compare an experiment sample against a control sample (CRISPRessoCompare).

    Only generated when the 'comparison' column is present in samples.csv and
    contains both 'experiment' and 'control' values. Every experiment sample is
    compared to every control sample (Cartesian product).
    """
    input:
        folder_1 = "crispresso/CRISPResso_on_{exp_name}",
        folder_2 = "crispresso/CRISPResso_on_{ctrl_name}",
    output:
        directory("crispresso/CRISPRessoCompare_{exp_name}_vs_{ctrl_name}"),
    params:
        sample_1_name = lambda wc: wc.exp_name,
        sample_2_name = lambda wc: wc.ctrl_name,
    log:
        "logs/{exp_name}_vs_{ctrl_name}.compare.log",
    wrapper:
        wrapper_path("CRISPR/crispresso-compare")


rule crispresso_aggregate:
    """Aggregate all CRISPResso runs into a single combined report.

    Collects every CRISPResso_on_{sample_name} directory produced by the
    crispresso rule and passes them to CRISPRessoAggregate, which produces a
    unified HTML report and summary plots across all samples.
    """
    input:
        crispresso_dirs = get_all_crispresso_outputs,
    output:
        directory("crispresso/CRISPRessoAggregate_on_all"),
    params:
        name = "all",
    threads: 4
    log:
        "logs/aggregate.log",
    wrapper:
        wrapper_path("CRISPR/crispresso-aggregate")
