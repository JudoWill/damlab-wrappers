# Deletion-pattern split + per-deletion samtools consensus rules
from os.path import join


def get_deletion_bam(wildcards):
    return f'consensus_split/{wildcards.sample}/{wildcards.category}.sorted.bam'


rule prepare_deletion_split_csv:
    """Transform deletion_reads.csv into a split-compatible query_name/category CSV.

    Filters:
      1. min_reference_size: drop reads whose reference span is below this threshold.
      2. min_reads_per_deletion: drop categories with fewer than this many reads.
    """
    input:
        reads='deletion_detection/{sample}.deletion_reads.csv'
    output:
        csv='deletion_split/{sample}.split_labels.csv'
    params:
        min_reference_size=config.get('MIN_REFERENCE_SIZE', 0),
        min_reads_per_deletion=config.get('MIN_READS_PER_DELETION', 10)
    run:
        import pandas as pd

        df = pd.read_csv(input.reads)

        # Filter 1: minimum reference span
        min_ref = params.min_reference_size
        if min_ref > 0:
            df = df[df['reference_end'] - df['reference_start'] >= min_ref]

        # Build category: empty deletions → "intact"; sanitize semicolons
        df['category'] = (
            df['deletions']
            .fillna('')
            .str.strip()
            .replace('', 'intact')
            .str.replace(';', '__', regex=False)
        )
        # Rows with an empty string after strip that weren't caught by replace
        df.loc[df['category'] == '', 'category'] = 'intact'

        # Filter 2: minimum reads per deletion category
        min_reads = params.min_reads_per_deletion
        if min_reads > 0:
            counts = df['category'].value_counts()
            keep = counts[counts >= min_reads].index
            df = df[df['category'].isin(keep)]

        # Write split-compatible CSV
        df[['read_name', 'category']].rename(
            columns={'read_name': 'query_name'}
        ).to_csv(output.csv, index=False)


rule namesort_for_deletion_split:
    """Name-sort the aligned BAM so the split wrapper groups supplementary alignments."""
    input:
        'aligned/{sample}.sorted.bam'
    output:
        temp('deletion_split/{sample}.namesorted.bam')
    log:
        'deletion_split/{sample}.namesort.log'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/samtools/collate"


checkpoint split_by_deletion_pattern:
    """Split name-sorted BAM into per-deletion-pattern BAM files."""
    input:
        bam='deletion_split/{sample}.namesorted.bam',
        csv_file='deletion_split/{sample}.split_labels.csv'
    output:
        directory('deletion_split/{sample}/')
    params:
        split_by='csv',
        output_format='bam',
        sample_name=lambda wildcards: wildcards.sample
    log:
        'deletion_split/{sample}.split.log'
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/cigarmath/split/"


rule sort_deletion_bam:
    """Coordinate-sort each per-deletion BAM for samtools consensus."""
    input:
        'deletion_split/{sample}/{category}.bam'
    output:
        temp('consensus_split/{sample}/{category}.sorted.bam')
    log:
        'consensus_split/{sample}/{category}.sort.log'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/samtools/sort"


def get_all_consensus_outputs(wildcards):
    """Return all samtools consensus outputs for a sample after the checkpoint.

    The 'unclassified' category is excluded: it contains reads that were not
    in the split CSV (non-NFL short/partial reads).
    """
    checkpoint_output = checkpoints.split_by_deletion_pattern.get(
        sample=wildcards.sample
    ).output[0]
    categories = [
        c for c in glob_wildcards(
            join(checkpoint_output, "{category}.bam")
        ).category
        if c != 'unclassified'
    ]

    return expand(
        'consensus_split/{sample}/{category}.consensus.fa',
        sample=wildcards.sample,
        category=categories
    )


rule samtools_consensus_per_deletion:
    """Run samtools consensus on each deletion-pattern BAM."""
    input:
        get_deletion_bam
    output:
        consensus='consensus_split/{sample}/{category}.consensus.fa'
    params:
        mode='bayesian'
    threads: 4
    log:
        'consensus_split/{sample}/{category}.consensus.log'
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/samtools/consensus/"


rule concatenate_haplotypes_for_msa:
    """Concatenate the reference and all per-deletion consensus sequences into a single FASTA."""
    input:
        haplotypes=get_all_consensus_outputs,
        reference=get_reference_index
    output:
        'consensus_split/{sample}.pre_msa.fasta'
    log:
        'consensus_split/{sample}.pre_msa.log'
    shell:
        'cat {input.reference} {input.haplotypes} > {output} 2> {log}'


rule muscle_deletion_msa:
    """Align reference + all consensus sequences with MUSCLE."""
    input:
        'consensus_split/{sample}.pre_msa.fasta'
    output:
        'consensus_split/{sample}.msa.fasta'
    threads: 4
    log:
        'consensus_split/{sample}.msa.log'
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/MSA/muscle/"
