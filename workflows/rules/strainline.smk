# Strainline haplotype reconstruction rules
from os.path import join

rule bam_to_fastq:
    """Convert aligned BAM to FASTQ for strainline input"""
    input:
        'aligned/{sample}.sorted.bam'
    output:
        temp('strainline/{sample}.fastq')
    params:
        extra="-F 2308" # Skip any unmapped, not primary, and secondary alignments
    log:
        'strainline/{sample}.bam2fastq.log'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/samtools/fastq/interleaved"

rule fastq_to_fasta:
    """Convert aligned BAM to FASTQ for strainline input"""
    input:
        'strainline/{sample}.fastq'
    output:
        temp('strainline/{sample}.fasta')
    params:
        command="seq"
    log:
        'strainline/{sample}.fastq2fasta.log'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/seqkit"

rule strainline:
    input:
        'strainline/{sample}.fasta'
    output:
        haplotypes='strainline/{sample}.haplotypes.fa'
    params:
        prefix=config.get('STRAINLINE_PREFIX', join(WORKFLOW_DIR, '../strainline/venv')),
        platform="ont"
    threads: 4
    log:
        'strainline/{sample}.strainline.log'
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/strainline/strainline/"
