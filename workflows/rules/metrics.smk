rule metrics:
    input:
        reads = 'aligned/{sample}.mm2.bam'
    output:
        'metrics/{sample}.hivmetrics.yaml'
    params:
        sample_name = lambda wildcards: wildcards.sample
    log:
        'metrics/{sample}.hivmetrics.log'
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/multiqc/hivmetrics/"

rule depth:
    input:
        reads = 'aligned/{sample}.mm2.bam'
    output:
        'metrics/{sample}.depth.txt'
    log:
        'metrics/{sample}.depth.log'
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/cigarmath/pileup/"

rule no_demux:
    output:
        'metrics/{sample}.nodemux.txt'
    params:
        sample_name = lambda wildcards: wildcards.sample,
        demux = 'No Reads'
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/multiqc/generic/"

