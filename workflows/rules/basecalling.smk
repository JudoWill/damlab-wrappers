from os.path import join
from pathlib import Path
import glob

def get_pod5_files(wildcards):
    """Get pod5 files from config or auto-detect"""
    # First check if pod5_path is specified in config
    pod5_path = config.get('pod5_path')
    if pod5_path:
        if Path(pod5_path).exists():
            return pod5_path
        else:
            raise ValueError(f'pod5_path specified in config ({pod5_path}) does not exist')
    
    # Fallback to original auto-detection logic
    if Path('converted/converted.pod5').exists():
        return 'converted/converted.pod5'
    elif Path('archive_pod5/').exists():
        return 'archive_pod5/'
    elif Path('pod5/').exists():
        return 'pod5/'
    
    raise ValueError('No pod5 files found. Specify pod5_path in config or ensure pod5 files exist in converted/, archive_pod5/, or pod5/')

def get_basecalled_bam():
    """Get path to basecalled bam based on DORADO_MODE"""
    if config.get('DORADO_MODE', 'duplex') == 'duplex':
        return 'duplex/basecalled.bam'
    else:
        return 'simplex/basecalled.bam'

# Helper functions for scatter mode
def use_scatter_mode():
    """Check if scatter mode is enabled"""
    return config.get('DORADO_SCATTER', False)

def get_pod5_files_list(wildcards):
    """Get list of individual pod5 files for scatter mode"""
    pod5_path = config.get('pod5_path')
    if pod5_path:
        path_obj = Path(pod5_path)
        if path_obj.exists():
            if path_obj.is_dir():
                # Find all pod5 files in the directory
                pod5_files = list(path_obj.glob('*.pod5'))
                return [str(f) for f in pod5_files]
            else:
                # Single file
                return [str(path_obj)]
    
    # Fallback to auto-detection
    if Path('converted/converted.pod5').exists():
        return ['converted/converted.pod5']
    elif Path('archive_pod5/').exists():
        pod5_files = list(Path('archive_pod5/').glob('*.pod5'))
        return [str(f) for f in pod5_files]
    elif Path('pod5/').exists():
        pod5_files = list(Path('pod5/').glob('*.pod5'))
        return [str(f) for f in pod5_files]
    
    raise ValueError('No pod5 files found for scatter mode')

def get_pod5_path_from_basename(wildcards):
    """Get full pod5 path from basename wildcard"""
    # Read the checkpoint to get the mapping
    checkpoint_output = checkpoints.pod5_discovery.get(**wildcards).output[0]
    
    with open(checkpoint_output, 'r') as f:
        pod5_files = [line.strip() for line in f if line.strip()]
    
    # Find the file matching the basename
    for pod5_file in pod5_files:
        if Path(pod5_file).stem == wildcards.basename:
            return pod5_file
    
    raise ValueError(f"Could not find pod5 file with basename {wildcards.basename}")

def aggregate_scattered_bams(wildcards):
    """Aggregate scattered BAM files after checkpoint resolution"""
    # Get the checkpoint output
    checkpoint_output = checkpoints.pod5_discovery.get(**wildcards).output[0]
    
    # Read the list of pod5 files from the checkpoint
    with open(checkpoint_output, 'r') as f:
        pod5_files = [line.strip() for line in f if line.strip()]
    
    # Determine mode (duplex or simplex)
    mode = config.get('DORADO_MODE', 'duplex')
    
    # Generate scattered BAM paths
    scattered_bams = []
    for pod5_file in pod5_files:
        # Extract basename without extension
        basename = Path(pod5_file).stem
        scattered_bams.append(f'{mode}/scattered/{basename}.bam')
    
    return scattered_bams

# Checkpoint for scatter mode: discover pod5 files
checkpoint pod5_discovery:
    output:
        'pod5_discovery/pod5_files.txt'
    run:
        # Get list of pod5 files
        pod5_files = get_pod5_files_list(wildcards)
        
        # Create output directory
        Path('pod5_discovery').mkdir(parents=True, exist_ok=True)
        
        # Write list to file
        with open(output[0], 'w') as f:
            for pod5_file in pod5_files:
                f.write(f'{pod5_file}\n')

# Scattered basecalling rules (one job per pod5 file)
rule pod5_duplexing_scattered:
    input:
        pod = get_pod5_path_from_basename
    output:
        temp('duplex/scattered/{basename}.bam')
    params:
        model = config.get('DORADO_MODEL', config.get('DUPLEX_MODEL', 'dna_r10.4.1_e8.2_5khz_stereo@v1.3')),
        models_directory = config.get('MODEL_ROOT'),
    threads: config.get('SCATTER_THREADS', workflow.cores)
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/dorado/duplex/"

rule pod5_simplex_scattered:
    input:
        pod = get_pod5_path_from_basename
    output:
        temp('simplex/scattered/{basename}.bam')
    params:
        model = config.get('DORADO_MODEL', config.get('SIMPLEX_MODEL', 'dna_r9.4.1_e8_sup@v3.6')),
        models_directory = config.get('MODEL_ROOT'),
    threads: config.get('SCATTER_THREADS', workflow.cores)
    wrapper:
        "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/dorado/simplex/"

# Merge scattered BAMs (only runs in scatter mode)
if use_scatter_mode():
    rule merge_duplex_scattered:
        input:
            aggregate_scattered_bams
        output:
            protected('duplex/basecalled.bam')
        log:
            'duplex/merge.log'
        threads: config.get('MERGE_THREADS', 16)
        wrapper:
            f"{config.get('snakemake_wrapper_tag', 'v8.1.1')}/bio/samtools/merge"

    rule merge_simplex_scattered:
        input:
            aggregate_scattered_bams
        output:
            protected('simplex/basecalled.bam')
        log:
            'simplex/merge.log'
        threads: config.get('MERGE_THREADS', 16)
        wrapper:
            f"{config.get('snakemake_wrapper_tag', 'v8.1.1')}/bio/samtools/merge"

# Non-scattered basecalling rules (original behavior, only runs when scatter mode is disabled)
if not use_scatter_mode():
    rule pod5_duplexing:
        input:
            pod = get_pod5_files
        output:
            protected('duplex/basecalled.bam')
        params:
            model = config.get('DORADO_MODEL', config.get('DUPLEX_MODEL', 'dna_r10.4.1_e8.2_5khz_stereo@v1.3')),
            models_directory = config.get('MODEL_ROOT'),
        threads: workflow.cores
        wrapper:
            "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/dorado/duplex/"

    rule pod5_simplex:
        input:
            pod = get_pod5_files
        output:
            protected('simplex/basecalled.bam')
        params:
            model = config.get('DORADO_MODEL', config.get('SIMPLEX_MODEL', 'dna_r9.4.1_e8_sup@v3.6')),
            models_directory = config.get('MODEL_ROOT'),
        threads: workflow.cores
        wrapper:
            "https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/dorado/simplex/"

