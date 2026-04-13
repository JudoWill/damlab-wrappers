"""proviral_crispr.smk — CRISPResso automation pipeline.

Automates CRISPResso (and optional CRISPRessoCompare) analysis from flexible
input types: paired or single-end FASTQ, or a BAM file (with optional region
slicing via cigarmath/slice).

Config keys (config.yaml or run.meta.yaml):
    samples_csv     : path to samples CSV (default: samples.csv)
    MIN_DELETION_SIZE       : optional; min deletion length for deletion_block_detection (default: 50)
    DELETION_MERGE_DISTANCE : optional; merge nearby deletion blocks (default: 10)
    DEBUG_DELETION_QUERY    : optional; if true (default), wrapper logs query param resolution to the rule log
    CRISPRESSO_EXTRA        : optional; extra CLI args for CRISPResso (crispresso rule). Default: --disable_guardrails
    MIN_READS_FOR_CRISPRESSO : optional; BAM-derived FASTQs (bam_to_fastq / slice) need at least this many reads (line_count // 4) for CRISPResso; default 10. CSV fastq_r1 mode is not gated.
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
                            cigarmath/deletion_block_detection as params.deletion_query. BAM samples only.
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
    "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main"
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


def _normalize_sample_table_columns(df):
    """Strip UTF-8 BOM and surrounding whitespace from CSV headers (common Excel/export issue)."""
    out = df.copy()
    out.columns = (
        out.columns.astype(str)
        .str.replace("\ufeff", "", regex=False)
        .str.strip()
    )
    return out


SAMPLES = _normalize_sample_table_columns(
    pd.read_csv(config.get("samples_csv", "samples.csv"), encoding="utf-8-sig")
)

# Align sample_name with Snakemake wildcards (trim accidental spaces from spreadsheet exports)
if "sample_name" in SAMPLES.columns:
    SAMPLES["sample_name"] = SAMPLES["sample_name"].astype(str).str.strip()

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


READ_GATE_NO_BAM_DERIVED = "qc/.read_gate_no_bam_derived"


def bam_derived_fastq_paths():
    """Paths to FASTQs this workflow generates from BAM (slice or bam2fastx), in CSV row order."""
    paths = []
    for _, row in SAMPLES.iterrows():
        mode = get_input_mode(row)
        if mode == "bam_region":
            paths.append(f"fastq/{row['sample_name']}.slice.fastq")
        elif mode == "bam":
            paths.append(f"fastq/{row['sample_name']}.bam.fastq")
    return paths


def read_gate_checkpoint_fastqs(wildcards):
    """Checkpoint inputs: all BAM-derived FASTQs, or a sentinel when there are none."""
    p = bam_derived_fastq_paths()
    return p if p else [READ_GATE_NO_BAM_DERIVED]


def _sample_name_from_bam_derived_fq(fq_path):
    base = os.path.basename(str(fq_path))
    if base.endswith(".slice.fastq"):
        return base[: -len(".slice.fastq")]
    if base.endswith(".bam.fastq"):
        return base[: -len(".bam.fastq")]
    raise ValueError(f"Unexpected BAM-derived FASTQ path: {fq_path}")


def _line_count_fast(path):
    n = 0
    with open(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n


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


def _cell_by_logical_column(row, logical_name: str):
    """Read one cell; match column by case-insensitive stripped name (handles odd CSV headers)."""
    if not hasattr(row, "index"):
        return None
    want = logical_name.strip().lower()
    for col in row.index:
        key = str(col).replace("\ufeff", "").strip().lower()
        if key == want:
            return row[col]
    return None


def _opt_deletion_query(wildcards):
    """``deletion_query`` cell as a trimmed string (always ``str()`` for pandas/numpy scalars)."""
    row = get_sample(wildcards)
    v = _cell_by_logical_column(row, "deletion_query")
    if v is None and hasattr(row, "get"):
        v = row.get("deletion_query")
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

def crispresso_eligible_samples(wildcards):
    """Samples that receive CRISPResso outputs: all fastq-mode rows plus BAM-derived rows that pass read gate."""
    ckpt = checkpoints.read_gate_aggregate.get()
    passed_path = ckpt.output.passed
    with open(passed_path) as fh:
        passed_bam = {ln.strip() for ln in fh if ln.strip()}
    eligible = []
    for _, row in SAMPLES.iterrows():
        name = row["sample_name"]
        mode = get_input_mode(row)
        if mode == "fastq" or name in passed_bam:
            eligible.append(name)
    return eligible


def get_all_crispresso_outputs(wildcards):
    return [
        f"crispresso/CRISPResso_on_{s}"
        for s in crispresso_eligible_samples(wildcards)
    ]


def get_all_compare_outputs(wildcards):
    if not HAS_COMPARISON:
        return []
    eligible = set(crispresso_eligible_samples(wildcards))
    return [
        f"crispresso/CRISPRessoCompare_{exp}_vs_{ctrl}"
        for exp, ctrl in itertools.product(EXPERIMENT_SAMPLES, CONTROL_SAMPLES)
        if exp in eligible and ctrl in eligible
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


rule read_gate_no_bam_derived:
    """Placeholder input when samples.csv has no BAM-derived FASTQs (checkpoint needs non-empty input)."""
    output:
        touch(READ_GATE_NO_BAM_DERIVED),


checkpoint read_gate_aggregate:
    """After BAM-derived FASTQs exist, list samples with enough reads for CRISPResso."""
    input:
        fastqs=read_gate_checkpoint_fastqs,
    output:
        passed="qc/passed_bam_samples.txt",
    log:
        "logs/read_gate_aggregate.log",
    run:
        min_reads = int(config.get("MIN_READS_FOR_CRISPRESSO", 10))
        out_path = Path(output.passed)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        raw = input.fastqs
        if isinstance(raw, str):
            paths = [raw]
        else:
            paths = list(raw)
        lines_out = []
        with open(log[0], "w") as logf:
            if (
                len(paths) == 1
                and str(paths[0]).endswith(".read_gate_no_bam_derived")
            ):
                logf.write("No BAM-derived FASTQs in samples.csv; gate passes none.\n")
            else:
                for fq in paths:
                    n_lines = _line_count_fast(fq)
                    n_reads = n_lines // 4
                    sample = _sample_name_from_bam_derived_fq(fq)
                    logf.write(
                        f"{sample}\t{fq}\tlines={n_lines}\treads={n_reads}\t"
                        f"min_reads={min_reads}\tpass={n_reads >= min_reads}\n"
                    )
                    if n_reads >= min_reads:
                        lines_out.append(sample)
        out_path.write_text("\n".join(lines_out) + ("\n" if lines_out else ""))


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
        debug_deletion_query=config.get("DEBUG_DELETION_QUERY", True),
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
        extra                        = config.get("CRISPRESSO_EXTRA", "--disable_guardrails"),
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
