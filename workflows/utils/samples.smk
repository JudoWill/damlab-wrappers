import csv
import yaml

def get_all_samples():
    """Get all unique sample names from SAMPLES DataFrame"""
    return list(SAMPLES['sample_name'].unique())

def get_sample_row(sample):
    """Get the row from SAMPLES DataFrame matching the sample name"""
    for _, row in SAMPLES.iterrows():
        if row['sample_name'] == sample:
            return row
    raise KeyError(f'Could not find row to match {sample}.')

def get_barcode_to_sample_map(wildcards):
    """Create mapping from barcode to output index for demux wrapper"""
    barcode_map = {}
    for idx, sample in enumerate(get_all_samples()):
        row = get_sample_row(sample)
        barcode_map[row['barcode']] = idx
    return barcode_map

def get_barcode_to_output(wildcards, input):
    """Create mapping from barcode to sample name for checkpoint demux"""
    with open(input['samples'], 'r') as f:
        samples = csv.DictReader(f)
        return {
            sample['barcode']: sample['sample_name'] for sample in samples
        }

def get_kit_name(wildcards, input):
    """Get kit name from runconfig, with SQK-NBD110-24 → SQK-NBD114-24 mapping"""
    with open(input['runconfig'], 'r') as f:
        runconfig = yaml.safe_load(f)
        kit = runconfig['library_prep_kits'][0]
        if kit == "SQK-NBD110-24":
            return "SQK-NBD114-24"
        return kit

