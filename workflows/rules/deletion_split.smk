# Deletion-pattern split + per-deletion Strainline rules
from os.path import join

HAS_DOWNSAMPLE = 'FILTLONG_TARGET_BASES' in config


def get_deletion_fasta(wildcards):
    if HAS_DOWNSAMPLE:
        return f'deletion_split/{wildcards.sample}/{wildcards.category}.downsampled.fa'
    return f'deletion_split/{wildcards.sample}/{wildcards.category}.fa'


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


if HAS_DOWNSAMPLE:
    _filtlong_extra = (
        f"--assembly {config['FILTLONG_ASSEMBLY']}"
        if config.get('FILTLONG_ASSEMBLY') else ""
    )

    rule downsample_deletion_split:
        """Downsample each per-deletion FASTA with filtlong before Strainline."""
        input:
            reads='deletion_split/{sample}/{category}.fa'
        output:
            temp('deletion_split/{sample}/{category}.downsampled.fastq')
        params:
            target_bases=config['FILTLONG_TARGET_BASES'],
            extra=_filtlong_extra
        log:
            'deletion_split/{sample}/{category}.filtlong.log'
        wrapper:
            f"{SNAKEMAKE_WRAPPER_TAG}/bio/filtlong"

    rule fq2fa_deletion_split:
        """Convert downsampled FASTQ back to FASTA for Strainline."""
        input:
            'deletion_split/{sample}/{category}.downsampled.fastq'
        output:
            'deletion_split/{sample}/{category}.downsampled.fa'
        params:
            command='fq2fa'
        log:
            'deletion_split/{sample}/{category}.fq2fa.log'
        wrapper:
            f"{SNAKEMAKE_WRAPPER_TAG}/bio/seqkit"


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
            join(checkpoint_output, "{category}.fa")
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
        get_deletion_fasta
    output:
        haplotypes='strainline_split/{sample}/{category}.haplotypes.fa'
    params:
        prefix=config.get('STRAINLINE_PREFIX', join(WORKFLOW_DIR, '../strainline/venv')),
        extra_params = '--minTrimmedLen 250 --minOvlpLen 250 --minSeedLen 500',
        platform='ont'
    threads: 4
    log:
        'strainline_split/{sample}/{category}.strainline.log'
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/strainline/strainline/"


rule concatenate_haplotypes_for_msa:
    """Concatenate the reference and all per-deletion haplotypes into a single FASTA."""
    input:
        haplotypes=get_all_strainline_split_outputs,
        reference=get_reference_index
    output:
        'strainline_split/{sample}.pre_msa.fasta'
    log:
        'strainline_split/{sample}.pre_msa.log'
    shell:
        'cat {input.reference} {input.haplotypes} > {output} 2> {log}'


rule muscle_deletion_msa:
    """Align reference + all haplotypes with MUSCLE."""
    input:
        'strainline_split/{sample}.pre_msa.fasta'
    output:
        'strainline_split/{sample}.msa.fasta'
    threads: 4
    log:
        'strainline_split/{sample}.msa.log'
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/MSA/muscle/"
