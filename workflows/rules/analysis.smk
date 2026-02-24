# Analysis rules for strainline haplotype reconstruction and deletion block detection
from os.path import join

rule bam_to_fastq:
    """Convert aligned BAM to FASTQ for strainline input"""
    input:
        bam='aligned/{sample}.sorted.bam'
    output:
        fastq=temp('analysis/{sample}.fastq')
    log:
        'analysis/{sample}.bam2fastq.log'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/samtools/fastq"

rule strainline:
    input:
        'analysis/{sample}.fastq'
    output:
        haplotypes='analysis/{sample}.haplotypes.fa'
    params:
        prefix=config.get('STRAINLINE_PREFIX', join(WORKFLOW_DIR, '../strainline/venv')),
        platform="ont"
    threads: 4
    log:
        'analysis/{sample}.strainline.log'
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/strainline/strainline/"

rule deletion_block_detection:
    input:
        bams='aligned/{sample}.sorted.bam'
    output:
        reads='analysis/{sample}.deletion_reads.csv',
        deletions='analysis/{sample}.deletion_blocks.csv',
        summary='analysis/{sample}.deletion_summary.yaml'
    params:
        min_deletion_size=config.get('MIN_DELETION_SIZE', 50),
        sample_name=lambda wildcards: wildcards.sample
    log:
        'analysis/{sample}.deletion_detection.log'
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/cigarmath/deletion_block_detection/"
