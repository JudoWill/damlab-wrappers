rule samtools_stats:
    input:
        bam = 'aligned/{sample}.sorted.bam',
    output:
        'samtools_stats/{sample}.txt',
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/samtools/stats"



rule run_report:
    input:
        expand('metrics/{sample_name}.hivmetrics.yaml', sample_name=SAMPLES['sample_name']),
        expand('metrics/{sample_name}.depth.txt', sample_name=SAMPLES['sample_name']),
        expand('samtools_stats/{sample_name}.txt', sample_name=SAMPLES['sample_name']),
        expand('strainline/{sample_name}.haplotypes.fa', sample_name=SAMPLES['sample_name']),
        expand('deletion_detection/{sample_name}.deletion_summary.yaml', sample_name=SAMPLES['sample_name'])
    output:
        report=f'qc/multiqc.html'
    params:
        prefix = f'{WORKFLOW_DIR}/../multiqc/venv',
        extra_args = '-v -f'
    log:
        f'qc/report.log'
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/multiqc/multiqc/"

