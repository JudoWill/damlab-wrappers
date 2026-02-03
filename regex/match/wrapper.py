"""Wrapper for regex pattern matching in FASTA/FASTQ files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import regex
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import yaml
import gzip
import os

if "snakemake" not in locals():
    import snakemake # type: ignore

def get_regex_patterns(params):
    """Get regex patterns from parameters"""
    if 'patterns' not in params:
        raise ValueError("No regex patterns provided in params. Please provide patterns as a dictionary in the 'patterns' parameter.")
    
    patterns = params['patterns']
    if not isinstance(patterns, dict):
        raise ValueError("Patterns must be provided as a dictionary where keys are column names and values are regex patterns")
    
    return {name: regex.compile(pattern, flags=regex.BESTMATCH+regex.IGNORECASE) 
            for name, pattern in patterns.items()}

def extract_matches(reg, text, both_strands=False):
    """Extract first match from text using regex pattern, optionally checking both strands"""
    # Try forward strand first
    match = reg.findall(text)
    if match:
        return match[0]
    
    # If no match and both_strands is True, try reverse complement
    if both_strands:
        rev_comp = str(Seq(text).reverse_complement())
        match = reg.findall(rev_comp)
        if match:
            return match[0]
    
    return None

def get_sequence_reader(file_path):
    """Determine the appropriate sequence reader based on file extension"""
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()
    
    # Handle gzipped files
    if ext == '.gz':
        _, base_ext = os.path.splitext(file_path[:-3])  # Remove .gz and get base extension
        base_ext = base_ext.lower()
        if base_ext in ['.fasta', '.fa', '.fna']:
            return lambda f: SeqIO.parse(gzip.open(f, 'rt'), 'fasta')
        elif base_ext in ['.fastq', '.fq']:
            return lambda f: SeqIO.parse(gzip.open(f, 'rt'), 'fastq')
        else:
            raise ValueError(f"Unsupported gzipped file format: {base_ext}")
    
    # Handle uncompressed files
    if ext in ['.fasta', '.fa', '.fna']:
        return lambda f: SeqIO.parse(f, 'fasta')
    elif ext in ['.fastq', '.fq']:
        return lambda f: SeqIO.parse(f, 'fastq')
    else:
        raise ValueError(f"Unsupported file format: {ext}")

def iterate_reads(file_path):
    """Iterate over sequence records in a file, handling compression and format automatically"""
    reader = get_sequence_reader(file_path)
    
    # Open file with appropriate mode based on compression
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as handle:
            yield from reader(file_path)
    else:
        with open(file_path, 'r') as handle:
            yield from reader(file_path)

def process_sequences(file_path, patterns, both_strands=False):
    """Process sequence file and extract matches for each pattern"""
    results = []
    
    for record in iterate_reads(file_path):
        row = {'read_name': record.id}
        
        # Extract matches for each pattern
        for name, pattern in patterns.items():
            match = extract_matches(pattern, str(record.seq), both_strands)
            row[name] = match
        
        results.append(row)
    
    return pd.DataFrame(results)

# Get input/output paths
seq_in_path = str(snakemake.input[0])
csv_out_path = str(snakemake.output[0])
metrics_file = snakemake.output.get('metrics', None)

# Get parameters
params = dict(snakemake.params)
both_strands = params.get('both_strands', False)
sample_name = params.get('sample_name', None)

# Get regex patterns
patterns = get_regex_patterns(params)

# Process sequence file
df = process_sequences(seq_in_path, patterns, both_strands)

# Write results to CSV
df.to_csv(csv_out_path, index=False)

# Generate and write metrics if requested
if metrics_file:
    metrics = {
        'total_reads': int(len(df)),  # Convert to native Python int
        'pattern_matches': {
            name: int(df[name].notna().sum())  # Convert to native Python int
            for name in patterns.keys()
        }
    }
    if sample_name:
        metrics['sample_name'] = sample_name
    
    with open(metrics_file, 'w') as f:
        f.write('# Regex pattern matching metrics\n')
        yaml.dump(metrics, f, default_flow_style=False) 