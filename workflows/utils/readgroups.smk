import csv
import yaml

def sample2readgroup_info(row):
    """Generate read group info from sample row and config.
    
    Uses config directly for run-level metadata (flow_cell_id, position, etc.)
    Row should contain sample-level info (sample_name, PCR_library, etc.)
    """
    # Per https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard

    info = {
        # SM Read-Group sample name
        'SM': row['sample_name'],
        
        # PL Read-Group platform (e.g. ILLUMINA, SOLID)
        'PL': 'ont',
        
        # PU Read-Group platform unit
        'PU': config.get('flow_cell_id', row.get('flow_cell_id')),
        
        # CN Read-Group sequencing center name
        'CN': config.get('position', row.get('position'))
        
    }
    
    #if row.get('description'):
    #    # DS Read-Group description
    #    info['DS'] = row['description']
    
    
    # LB Read-Group library
    #  Distinguishes reads from the same biological sample
    #  but from different PCRs
    info['LB'] = '.'.join([
        row['sample_name'],
        config.get('pcr_protocol', row.get('pcr_protocol')),
        str(row['PCR_library'])
    ])
    
    # ID Read-Group ID
    #  Distinguishes different runs of the same library
    #  ideally, unique in the world
    info['ID'] = '.'.join([
        str(info['LB']),
        config.get('run_id', row.get('run_id'))
    ])
    
    # PM Read-Group platform model
    #  Using this to indicate the basecalling-model.
    #  This also nicely denotes the chemistry.
    info['PM'] = config.get('SIMPLEX_MODEL', 'dna_r9.4.1_e8_sup@v3.6')
    if config.get('DUPLEX_MODEL'):
        info['PM'] += ',' + config['DUPLEX_MODEL']
    
    return info

def get_sample_rg_params(wildcards):
    """Get read group parameters for a sample from config and samples.csv"""
    # Read config (could be from run.meta.yaml or inline config)
    runconfig = config.copy()
    
    # Read samples.csv
    samples_csv = config.get('samples_csv', 'samples.csv')
    with open(samples_csv, 'r') as f:
        samples = csv.DictReader(f)
        for sample in samples:
            if sample['sample_name'] == wildcards.sample:
                # Merge sample info into runconfig
                runconfig.update(sample)
                return sample2readgroup_info(runconfig)
    
    raise ValueError(f"Sample {wildcards.sample} not found in samplesheet {samples_csv}")

def sample2readgroup_extra(wildcards):
    """Generate extra parameters string for picard AddOrReplaceReadGroups wrapper"""
    row = get_sample_row(wildcards['sample'])
    info = sample2readgroup_info(row)
    
    return ' '.join(f'-{arg} {val.replace("-","_")}' for arg, val in info.items())

