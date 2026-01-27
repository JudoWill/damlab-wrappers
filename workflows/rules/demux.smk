from os.path import join
import csv
import yaml

def get_barcode_to_output_from_config(wildcards, input):
    """Get barcode to output mapping from samples.csv"""
    mapping = {}
    for _, row in SAMPLES.iterrows():
        mapping[row['barcode']] = row['sample_name']
    return mapping

def get_kit_name_from_config(wildcards, input):
    """Get kit name from config, with SQK-NBD110-24 → SQK-NBD114-24 mapping"""
    kits = config.get('library_prep_kits', [])
    if isinstance(kits, list) and len(kits) > 0:
        kit = kits[0]
        if kit == "SQK-NBD110-24":
            return "SQK-NBD114-24"
        return kit
    raise ValueError("library_prep_kits not found in config or is empty")

# Checkpoint for demuxing from basecalled bam (pod5 input mode)
rule demux_run:
    input:
        reads = get_basecalled_bam(),
        
        # Not inputs, but needed for knowing when to re-run
        _runconfig = config.get('run_config', 'run.meta.yaml'),
        _samples = config.get('samples_csv', 'samples.csv')
    output:
        output_direc = directory('demux/'),
        output_files = temp(expand('demux/{sample}.bam', sample=get_all_samples()))
    params:
        kit_name = get_kit_name_from_config,
        barcode_to_output = get_barcode_to_output_from_config
    threads: int(workflow.cores * 0.75) + 1
    log:
        'demux.log'
    wrapper:
        #"https://raw.githubusercontent.com/DamLabResources/damlab-wrappers/refs/heads/main/dorado/demux/"
        "file:/home/jupyter-will/repos/damlab-wrappers/dorado/demux"
