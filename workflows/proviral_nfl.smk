from os.path import join
import pandas as pd
from pathlib import Path

# Support both run.meta.yaml and config.yaml for backwards compatibility
# Default to run.meta.yaml (original format), can override with --configfile config.yaml
configfile: 'run.meta.yaml'

# Constants with defaults for backwards compatibility
SNAKEMAKE_WRAPPER_TAG = config.get('snakemake_wrapper_tag', 'v8.1.1')
WORKFLOW_DIR = workflow.basedir
DL_PREFIX = config.get("damlab_prefix", join(WORKFLOW_DIR, "../../damlab-wrappers"))

# Load samples CSV
SAMPLES = pd.read_csv(config.get('samples_csv', 'samples.csv'))

# Include schema validation
include: join(WORKFLOW_DIR, "schemas", "validate.smk")

# Include utility modules
include: join(WORKFLOW_DIR, "utils/samples.smk")
include: join(WORKFLOW_DIR, "utils/readgroups.smk")

# Include rule modules
include: join(WORKFLOW_DIR, "rules/basecalling.smk")
include: join(WORKFLOW_DIR, "rules/demux.smk")
include: join(WORKFLOW_DIR, "rules/alignment.smk")
include: join(WORKFLOW_DIR, "rules/metrics.smk")
include: join(WORKFLOW_DIR, "rules/analysis.smk")
include: join(WORKFLOW_DIR, "rules/reporting.smk")

def get_final_bam_paths(wildcards):
    """Get all final sorted bam paths"""
    wanted = []
    for sample in get_all_samples():
        wanted.append(f'aligned/{sample}.sorted.bam')
    return wanted

def get_all_outputs(wildcards):
    """Get all final outputs for rule all"""
    outputs = get_final_bam_paths(wildcards)
    outputs.extend([
        'qc/multiqc.html',
        #'qc/multiqc_data.zip'
    ])
    # Add analysis outputs
    for sample in get_all_samples():
        outputs.append(f'analysis/{sample}.haplotypes.fa')
        outputs.append(f'analysis/{sample}.deletion_summary.yaml')
    return outputs

# Determine input mode
def get_input_mode():
    """Determine if we're using pod5 input or pre-demuxed bam input"""
    if config.get('demuxed_bam_path'):
        return 'demuxed'
    elif config.get('pod5_path') or Path('pod5/').exists() or Path('archive_pod5/').exists() or Path('converted/converted.pod5').exists():
        return 'pod5'
    else:
        # Default to pod5 mode if nothing specified
        return 'pod5'

rule all:
    input:
        get_all_outputs

rule all_pod5:
    """Entry point for pod5 input mode"""
    input:
        get_all_outputs

rule all_demuxed:
    """Entry point for pre-demuxed bam input mode"""
    input:
        get_all_outputs

