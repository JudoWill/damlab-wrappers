# Strainline haplotype reconstruction rules
from os.path import join

rule bam_to_fasta:
    """Convert aligned BAM to FASTA for strainline input"""
    input:
        'aligned/{sample}.sorted.bam'
    output:
        temp('strainline/{sample}.fasta')
    params:
        outputtype="fasta"
    log:
        'strainline/{sample}.bam2fasta.log'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/samtools/fastx"

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
