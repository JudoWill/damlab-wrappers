from os.path import join
from pathlib import Path

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

