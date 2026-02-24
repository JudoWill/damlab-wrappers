# Deletion block detection rules

rule deletion_block_detection:
    input:
        bams='aligned/{sample}.sorted.bam'
    output:
        reads='deletion_detection/{sample}.deletion_reads.csv',
        deletions='deletion_detection/{sample}.deletion_blocks.csv',
        summary='deletion_detection/{sample}.deletion_summary.yaml'
    params:
        min_deletion_size=config.get('MIN_DELETION_SIZE', 50),
        sample_name=lambda wildcards: wildcards.sample
    log:
        'deletion_detection/{sample}.deletion_detection.log'
    wrapper:
        "https://raw.githubusercontent.com/JudoWill/damlab-wrappers/refs/heads/main/cigarmath/deletion_block_detection/"
