"""proviral_reconstruction.smk — Proviral haplotype reconstruction pipeline.

Automates long-read haplotype reconstruction from nanopore sequencing data
using Strainline, with optional downsampling, VADR annotation, GenBank file
generation, and phylogenetic tree generation grouped by sample cohort.

Config keys (run.meta.yaml):
    samples_csv          : path to samples CSV (default: samples.csv)
    strainline_prefix    : (required) path to Strainline installation directory
    reference_fasta      : (required) reference FASTA for clipqs alignment/orientation
    downsample_size      : (optional) if set, filtlong target_bases for downsampling
    filtlong_reference   : (optional) FASTA reference for filtlong --assembly scoring

    VADR annotation — provide EXACTLY ONE of the three modes below:

    Mode A — pre-built model directory (no build step):
    annotation_genbank   : path to an existing vadr-build output directory.
                           The build and normalisation steps are skipped entirely.

    Mode B — fetch from NCBI (requires network, accession only):
    annotation_accession : GenBank/RefSeq accession (e.g. "K03455").
                           v-build.pl fetches the sequence and feature table from
                           NCBI at pipeline start.

    Mode C — build from local files (fully offline):
    annotation_accession : bare accession / model name (e.g. "K03455")
    annotation_local_genbank : path to a local GenBank flat file (.gb)
    annotation_local_fasta   : path to a local FASTA file (.fasta / .fa)
                           The vadr-genbank wrapper normalises the LOCUS name and
                           VERSION field before passing the files to vadr-build.

    If none of the above is set the annotation steps are skipped.

    Additional VADR options (Modes B and C only):
    annotation_group     : (optional) --group label for vadr-build (e.g. "HIV")
    annotation_subgroup  : (optional) --subgroup label (e.g. "HIV-1")

    GenBank output options (all modes):
    genbank_organism     : (optional) organism string for GenBank source records
                           (e.g. "Human immunodeficiency virus 1")

    damlab_prefix        : base location for damlab-wrappers. Can be:
                            - a local filesystem path  (e.g. /path/to/damlab-wrappers)
                            - a URL                    (e.g. https://raw.githubusercontent.com/...)
                           Default: https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main

samples.csv columns:
    sample_name  (required) unique run name
    path         (required) FASTQ or BAM file; BAM files are converted to FASTQ automatically
    tree_group   (optional) group label; samples sharing a label are aligned and
                            used to build a phylogenetic tree together
"""

import os
import itertools
import re
from pathlib import Path

import pandas as pd

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
# Validate required config
# ---------------------------------------------------------------------------

if "strainline_prefix" not in config:
    raise ValueError(
        "Config key 'strainline_prefix' is required. "
        "Set it in run.meta.yaml to the path of your Strainline installation."
    )

if "reference_fasta" not in config:
    raise ValueError(
        "Config key 'reference_fasta' is required. "
        "Set it in run.meta.yaml to a reference FASTA for clipqs alignment."
    )

STRAINLINE_PREFIX = config["strainline_prefix"]
REFERENCE_FASTA = config["reference_fasta"]

# ---------------------------------------------------------------------------
# Optional config flags
# ---------------------------------------------------------------------------

HAS_DOWNSAMPLE = "downsample_size" in config and config["downsample_size"]

# ---------------------------------------------------------------------------
# VADR annotation mode detection
#
#  Mode A — pre-built model dir   : annotation_genbank is set
#  Mode B — NCBI fetch            : annotation_accession set, no local files
#  Mode C — local files (offline) : annotation_accession + annotation_local_*
# ---------------------------------------------------------------------------

_ANNOT_GENBANK    = config.get("annotation_genbank", "")      # Mode A: pre-built dir
_ANNOT_ACCESSION  = config.get("annotation_accession", "")    # Modes B/C: accession
_ANNOT_LOCAL_GB   = config.get("annotation_local_genbank", "") # Mode C: local .gb
_ANNOT_LOCAL_FA   = config.get("annotation_local_fasta", "")   # Mode C: local .fasta

HAS_ANNOTATION     = bool(_ANNOT_GENBANK or _ANNOT_ACCESSION)
HAS_BUILD_STEP     = HAS_ANNOTATION and not _ANNOT_GENBANK
HAS_NORMALISE_STEP = HAS_BUILD_STEP and bool(_ANNOT_LOCAL_GB or _ANNOT_LOCAL_FA)

if HAS_NORMALISE_STEP and not _ANNOT_ACCESSION:
    raise ValueError(
        "annotation_accession must be set when using annotation_local_genbank "
        "or annotation_local_fasta (it provides the VADR model name)."
    )

if HAS_NORMALISE_STEP and not (_ANNOT_LOCAL_GB and _ANNOT_LOCAL_FA):
    raise ValueError(
        "Both annotation_local_genbank AND annotation_local_fasta must be "
        "provided together for local-file model building (Mode C)."
    )

# The VADR model directory used by vadr_annotate
if _ANNOT_GENBANK:
    VADR_MODEL_DIR = _ANNOT_GENBANK
elif _ANNOT_ACCESSION:
    VADR_MODEL_DIR = f"vadr_model/{_ANNOT_ACCESSION}"
else:
    VADR_MODEL_DIR = ""

# Normalised local-file staging directory (Mode C only)
_NORM_DIR = f"vadr_model/normalised/{_ANNOT_ACCESSION}" if HAS_NORMALISE_STEP else ""

# ---------------------------------------------------------------------------
# Load samples
# ---------------------------------------------------------------------------

SAMPLES = pd.read_csv(config.get("samples_csv", "samples.csv"))

# Validate required columns
for _col in ("sample_name", "path"):
    if _col not in SAMPLES.columns:
        raise ValueError(f"samples.csv must contain a '{_col}' column.")

HAS_TREE_GROUPS = (
    "tree_group" in SAMPLES.columns
    and SAMPLES["tree_group"].notna().any()
)

if HAS_TREE_GROUPS:
    TREE_GROUPS = SAMPLES["tree_group"].dropna().unique().tolist()
else:
    TREE_GROUPS = []

# Constrain sample_name to only match real sample names from the CSV.
# Without this, Snakemake can resolve e.g. fastq/HC69.downsampled.fastq
# via bam_to_fastq with sample_name=HC69.downsampled instead of using the
# downsample rule with sample_name=HC69.
wildcard_constraints:
    sample_name = "|".join(re.escape(s) for s in SAMPLES["sample_name"].tolist()),

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


def is_bam(path):
    """Return True if path points to a BAM file."""
    return str(path).lower().endswith(".bam")


def get_reads_input(wildcards):
    """Return the FASTQ path feeding into strainline (or filtlong).

    If the sample path is a BAM, the intermediate converted FASTQ is returned;
    otherwise the path from samples.csv is used directly.
    """
    row = get_sample(wildcards)
    if is_bam(row["path"]):
        return f"fastq/{wildcards.sample_name}.fastq"
    return row["path"]


def get_strainline_input(wildcards):
    """Return the path of reads to feed into Strainline.

    If downsampling is enabled, return the filtlong output; otherwise return
    the raw (or BAM-converted) FASTQ.
    """
    if HAS_DOWNSAMPLE:
        return f"fastq/{wildcards.sample_name}.downsampled.fastq"
    return get_reads_input(wildcards)


def get_bam2fastq_length_params(wildcards):
    """Return min/max query length params for bam_to_fastq from optional CSV columns."""
    row = get_sample(wildcards)
    params = {}
    if "min_query_length" in SAMPLES.columns and _notna(row.get("min_query_length")):
        params["min_query_length"] = int(row["min_query_length"])
    if "max_query_length" in SAMPLES.columns and _notna(row.get("max_query_length")):
        params["max_query_length"] = int(row["max_query_length"])
    return params


def get_samples_in_group(group):
    """Return sample names belonging to the given tree_group."""
    return SAMPLES[SAMPLES["tree_group"] == group]["sample_name"].tolist()


def get_clipped_fastas_for_group(wildcards):
    """Return all clipped FASTA paths for samples in a tree_group."""
    samples = get_samples_in_group(wildcards.tree_group)
    return [f"clipped/{s}.fa" for s in samples]


# ---------------------------------------------------------------------------
# Output collection helpers
# ---------------------------------------------------------------------------

def get_all_outputs(wildcards):
    all_outputs = []

    # Per-sample clipped haplotypes are always produced
    all_outputs += expand("clipped/{sample_name}.fa", sample_name=SAMPLES["sample_name"])

    # Optional per-sample VADR annotation + GenBank conversion
    if HAS_ANNOTATION:
        all_outputs += expand("vadr/{sample_name}", sample_name=SAMPLES["sample_name"])
        all_outputs += expand("gbk/{sample_name}.gbf", sample_name=SAMPLES["sample_name"])

    # Optional per-tree-group phylo outputs
    if HAS_TREE_GROUPS:
        all_outputs += expand("msa/{tree_group}.fa",      tree_group=TREE_GROUPS)
        all_outputs += expand("trees/{tree_group}.nwk",   tree_group=TREE_GROUPS)
        all_outputs += expand("figures/{tree_group}.png", tree_group=TREE_GROUPS)

    return all_outputs


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

rule all:
    input:
        get_all_outputs


# ---------------------------------------------------------------------------
# 1. BAM → FASTQ conversion (only when input is a BAM)
# ---------------------------------------------------------------------------

rule bam_to_fastq:
    """Convert a BAM file to FASTQ (cigarmath/bam2fastx).

    Triggered automatically when a sample's 'path' column points to a .bam file.
    Full read sequences are exported without region slicing.

    Optional samples.csv columns:
        min_query_length : discard reads shorter than this value (bp)
        max_query_length : discard reads longer than this value (bp)
    """
    input:
        lambda wc: get_sample(wc)["path"],
    output:
        fastq = "fastq/{sample_name}.fastq",
    params:
        mapped_only      = True,
        primary_only     = True,
        unique_only      = True,
        **get_bam2fastq_length_params,
    log:
        "logs/{sample_name}.bam2fastq.log",
    wrapper:
        wrapper_path("cigarmath/bam2fastx")


# ---------------------------------------------------------------------------
# 2. (Optional) Downsampling with filtlong
# ---------------------------------------------------------------------------

if HAS_DOWNSAMPLE:
    _filtlong_extra = ""
    if "filtlong_reference" in config and config["filtlong_reference"]:
        _filtlong_extra = f"--assembly {config['filtlong_reference']}"

    rule downsample:
        """Downsample reads with filtlong (snakemake-wrappers bio/filtlong).

        Uses 'downsample_size' from config as the target_bases limit.
        If 'filtlong_reference' is set in config, it is passed as --assembly
        so filtlong can up-weight reads that match the reference.
        """
        input:
            reads = get_reads_input,
        output:
            "fastq/{sample_name}.downsampled.fastq",
        params:
            target_bases = config["downsample_size"],
            extra        = _filtlong_extra,
        log:
            "logs/{sample_name}.filtlong.log",
        wrapper:
            "v7.6.0/bio/filtlong"


# ---------------------------------------------------------------------------
# 3. Strainline haplotype reconstruction
# ---------------------------------------------------------------------------

rule strainline_reconstruct:
    """Reconstruct haplotypes from long reads using Strainline.

    Produces a FASTA of assembled haplotype sequences per sample.
    Requires 'strainline_prefix' in config pointing to the Strainline
    installation directory.
    """
    input:
        get_strainline_input,
    output:
        haplotypes = "strainline/{sample_name}/haplotypes.fa",
    params:
        prefix   = STRAINLINE_PREFIX,
        platform = "ont",
    threads: 8
    log:
        "logs/{sample_name}.strainline.log",
    wrapper:
        wrapper_path("strainline/strainline")


# ---------------------------------------------------------------------------
# 4. Strainline clipping (orient + trim to reference)
# ---------------------------------------------------------------------------

rule strainline_clip:
    """Clip and orient Strainline haplotypes against the reference (strainline/clipqs).

    Aligns each haplotype to the reference FASTA, trims unaligned ends, and
    reverse-complements sequences on the minus strand.  Sequences with
    insufficient reference coverage are discarded.
    """
    input:
        sequences = "strainline/{sample_name}/haplotypes.fa",
        reference = REFERENCE_FASTA,
    output:
        "clipped/{sample_name}.fa",
    params:
        include_reference = False,
    log:
        "logs/{sample_name}.clipqs.log",
    wrapper:
        wrapper_path("strainline/clipqs")


# ---------------------------------------------------------------------------
# 5a. (Optional) VADR normalisation + model build
#
#   Mode A: annotation_genbank is a pre-built dir → both rules skipped.
#   Mode B: annotation_accession only → vadr_build fetches from NCBI.
#   Mode C: annotation_accession + local files → vadr_genbank normalises
#           the GenBank and FASTA files, then vadr_build builds offline.
# ---------------------------------------------------------------------------

if HAS_NORMALISE_STEP:
    rule vadr_genbank:
        """Normalise a local GenBank/FASTA pair for VADR compatibility (Mode C).

        VADR's GenBank parser requires the LOCUS name to equal the bare
        accession and the VERSION line to have no .N suffix.  Files downloaded
        from NCBI often violate both rules (e.g. K03455.gb has LOCUS HIVHXB2CG).
        This rule rewrites a clean copy; the originals are never modified.
        """
        input:
            genbank = _ANNOT_LOCAL_GB,
            fasta   = _ANNOT_LOCAL_FA,
        output:
            genbank = f"{_NORM_DIR}/{Path(_ANNOT_LOCAL_GB).name}",
            fasta   = f"{_NORM_DIR}/{Path(_ANNOT_LOCAL_FA).name}",
        params:
            accession = _ANNOT_ACCESSION,
        log:
            "logs/vadr_genbank.log",
        wrapper:
            wrapper_path("hiv/vadr-genbank")

    rule vadr_build:
        """Build a VADR model from normalised local files (Mode C — offline)."""
        input:
            genbank = f"{_NORM_DIR}/{Path(_ANNOT_LOCAL_GB).name}",
            fasta   = f"{_NORM_DIR}/{Path(_ANNOT_LOCAL_FA).name}",
        output:
            directory(VADR_MODEL_DIR),
        params:
            accession = _ANNOT_ACCESSION,
            group     = config.get("annotation_group", ""),
            subgroup  = config.get("annotation_subgroup", ""),
        log:
            "logs/vadr_build.log",
        wrapper:
            wrapper_path("hiv/vadr-build")

elif HAS_BUILD_STEP:
    rule vadr_build:
        """Build a VADR model by fetching an accession from GenBank (Mode B).

        Runs once per pipeline execution and stores the model in
        vadr_model/<accession>/.  Requires network access.
        """
        output:
            directory(VADR_MODEL_DIR),
        params:
            accession = _ANNOT_ACCESSION,
            group     = config.get("annotation_group", ""),
            subgroup  = config.get("annotation_subgroup", ""),
        log:
            "logs/vadr_build.log",
        wrapper:
            wrapper_path("hiv/vadr-build")


# ---------------------------------------------------------------------------
# 5b. (Optional) VADR annotation
# ---------------------------------------------------------------------------

if HAS_ANNOTATION:
    def _get_vadr_mdir(wildcards):
        """Return the VADR model directory for vadr_annotate.

        When a build step is required, the rule-output directory is returned
        so Snakemake can track the dependency.  When a pre-built directory is
        provided it is used directly (no build dependency).
        """
        return VADR_MODEL_DIR

    rule vadr_annotate:
        """Annotate clipped haplotypes using VADR (hiv/vadr-annotate).

        Uses the VADR model library in VADR_MODEL_DIR (either pre-built via
        annotation_genbank or freshly built from annotation_accession).
        Produces per-sample annotation directories with passing/failing
        classifications and feature tables.

        mode='hiv' activates --alt_pass for all recoverable alert codes so
        that proviral sequences with frameshifts, internal stops, etc. are
        still classified as passing.
        """
        input:
            sequences = "clipped/{sample_name}.fa",
            mdir      = _get_vadr_mdir,
        output:
            directory("vadr/{sample_name}"),
        params:
            mode         = "hiv",
            noseqnamemax = True,
        threads: 4
        log:
            "logs/{sample_name}.vadr_annotate.log",
        wrapper:
            wrapper_path("hiv/vadr-annotate")

    rule vadr_tbl2gbk:
        """Convert VADR passing sequences + feature table to GenBank format (hiv/vadr-tbl2gbk).

        Passes the vadr_annotate output directory directly; the wrapper locates
        the *.vadr.pass.fa and *.vadr.pass.tbl files within it automatically.
        This avoids a MissingInputException that would occur if specific files
        inside the directory() output were referenced directly.
        """
        input:
            vadr_dir = "vadr/{sample_name}",
        output:
            gbf = "gbk/{sample_name}.gbf",
            sqn = "gbk/{sample_name}.sqn",
        params:
            organism = config.get("genbank_organism", ""),
        log:
            "logs/{sample_name}.tbl2gbk.log",
        wrapper:
            wrapper_path("hiv/vadr-tbl2gbk")


# ---------------------------------------------------------------------------
# 6. (Optional) MSA per tree_group
# ---------------------------------------------------------------------------

if HAS_TREE_GROUPS:
    rule concat_group_fastas:
        """Concatenate clipped haplotypes for all samples in a tree_group.

        Produces a single multi-FASTA that is then aligned by MUSCLE.
        Each haplotype header is prefixed with the sample name for
        traceability in the tree.
        """
        input:
            get_clipped_fastas_for_group,
        output:
            "msa/{tree_group}.combined.fa",
        params:
            samples = lambda wc: get_samples_in_group(wc.tree_group),
        run:
            with open(output[0], "w") as out_fh:
                for fa_path, sample in zip(input, params.samples):
                    with open(fa_path) as in_fh:
                        for line in in_fh:
                            if line.startswith(">"):
                                out_fh.write(f">{sample}__{line[1:]}")
                            else:
                                out_fh.write(line)

    rule msa_by_group:
        """Multiple sequence alignment of haplotypes within a tree_group (MSA/muscle).

        Aligns the combined multi-FASTA for a group so that FastTree can infer
        a phylogeny.
        """
        input:
            "msa/{tree_group}.combined.fa",
        output:
            "msa/{tree_group}.fa",
        threads: 4
        log:
            "logs/{tree_group}.muscle.log",
        wrapper:
            wrapper_path("MSA/muscle")

    rule fasttree_by_group:
        """Infer a maximum-likelihood phylogeny per tree_group (phylo/FastTree).

        Uses the GTR+Gamma model, appropriate for HIV/retroviral sequences.
        """
        input:
            "msa/{tree_group}.fa",
        output:
            "trees/{tree_group}.nwk",
        params:
            gtr   = True,
            gamma = True,
            nt    = True,
        threads: 4
        log:
            "logs/{tree_group}.fasttree.log",
        wrapper:
            wrapper_path("phylo/FastTree")

    rule phytreeviz_by_group:
        """Visualize the phylogenetic tree per tree_group (phylo/phytreeviz).

        Produces a publication-quality PNG of the inferred tree.
        """
        input:
            "trees/{tree_group}.nwk",
        output:
            "figures/{tree_group}.png",
        params:
            format              = "newick",
            width               = 12,
            height              = 10,
            dpi                 = 300,
            show_branch_support = True,
        log:
            "logs/{tree_group}.phytreeviz.log",
        wrapper:
            wrapper_path("phylo/phytreeviz")
