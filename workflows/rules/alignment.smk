from os.path import join

def get_demuxed_reads(wildcards):
    """Get path to demuxed reads - supports both checkpoint and direct demux"""
    # Check if we're using pre-demuxed bam input mode
    if config.get('demuxed_bam_path'):
        # Using pre-demuxed input mode - use_demuxed_bam rule will provide the files
        return f'demux/{wildcards.sample}.bam'
    else:
        # Using pod5 input mode - check if checkpoint exists
        try:
            _ = checkpoints.demux_run.get().output
            return f'demux/{wildcards.sample}.bam'
        except:
            # Checkpoint not completed yet, will be resolved at runtime
            # Fall back to bam_demuxing rule output
            return f'demux/{wildcards.sample}.bam'

def get_reference_index(wildcards):
    """Get reference index path from config"""
    reference = config.get('REFERENCE')
    if not reference:
        raise ValueError("REFERENCE must be specified in config")
    return reference

rule add_read_groups:
    input:
        reads = get_demuxed_reads,
    output:
        temp('labeled/{sample}.rg.bam')
    params:
        extra = sample2readgroup_extra
    log:
        'labeled/{sample}.log'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/picard/addorreplacereadgroups"

rule align_reads:
    input:
        calls = 'labeled/{sample}.rg.bam',
        index = get_reference_index
    output:
        temp('aligned/{sample}.mm2.bam')
    log:
        'aligned/{sample}.log'
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/dorado/aligner/"

rule sort_reads:
    input:
        'aligned/{sample}.mm2.bam'
    output:
        'aligned/{sample}.sorted.bam'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/sambamba/sort"

rule index_reads:
    input:
        'aligned/{sample}.sorted.bam'
    output:
        'aligned/{sample}.sorted.bam.bai'
    wrapper:
        f"{SNAKEMAKE_WRAPPER_TAG}/bio/sambamba/index"

