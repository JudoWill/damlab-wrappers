# Deletion-pattern split + per-deletion Strainline rules
from os.path import join


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
    """Split name-sorted BAM into per-deletion-pattern FASTA files."""
    input:
        bam='deletion_split/{sample}.namesorted.bam',
        csv_file='deletion_split/{sample}.split_labels.csv'
    output:
        directory('deletion_split/{sample}/')
    params:
        split_by='csv',
        output_format='fasta',
        sample_name=lambda wildcards: wildcards.sample
    log:
        'deletion_split/{sample}.split.log'
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/cigarmath/split/"


def get_all_strainline_split_outputs(wildcards):
    """Return all strainline haplotype outputs for a sample after the checkpoint.

    The 'unclassified' category is excluded: it contains reads that were not
    in the split CSV (non-NFL short/partial reads) and will cause Strainline
    to hang or time out when the seed read is much shorter than a full-length
    genome.
    """
    checkpoint_output = checkpoints.split_by_deletion_pattern.get(
        sample=wildcards.sample
    ).output[0]
    categories = [
        c for c in glob_wildcards(
            join(checkpoint_output, "{category}.fasta")
        ).category
        if c != 'unclassified'
    ]
    return expand(
        'strainline_split/{sample}/{category}.haplotypes.fa',
        sample=wildcards.sample,
        category=categories
    )


rule strainline_per_deletion:
    """Run Strainline haplotype reconstruction on each deletion-pattern FASTA."""
    input:
        'deletion_split/{sample}/{category}.fasta'
    output:
        haplotypes='strainline_split/{sample}/{category}.haplotypes.fa'
    params:
        prefix=config.get('STRAINLINE_PREFIX', join(WORKFLOW_DIR, '../strainline/venv')),
        platform='ont'
    threads: 4
    log:
        'strainline_split/{sample}/{category}.strainline.log'
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/strainline/strainline/"


rule all_deletion_split_strainline:
    """Aggregate all per-deletion Strainline outputs for a sample."""
    input:
        get_all_strainline_split_outputs
    output:
        touch('strainline_split/{sample}.done')
