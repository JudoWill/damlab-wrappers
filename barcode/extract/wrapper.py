"""Wrapper for extracting barcodes and UMIs from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.1.0"

import pysam
import regex
from contextlib import contextmanager
import yaml

if "snakemake" not in locals():
    import snakemake # type: ignore


def _build_siv_barcode(strain, i, d, s):

    left = 'cctcctccccctccaggactagcataa'
    right = 'atggaagaaagacctccagaaaatga'
    # Taken from https://journals.asm.org/doi/full/10.1128/jvi.01420-19
    # Figure 5
    if strain == 'SIVmac239M':
        inner_left =  'ACGCGGGCTACC'
        inner_right = 'TTGCAAGCGCGT'
    elif strain == 'SIVmac239m2':
        inner_left =  'ACGCGAGCCATG'
        inner_right = 'CGATGCGCGCGT'
    elif strain == 'SIVmac239Opt5M':
        inner_left =  'ACGCGGACGGCA'
        inner_right = 'AGCGTGCCGCGT'
    elif strain == 'bSHIV1054M':
        inner_left =  'ACGCGGGCTACC'
        inner_right = 'GTGCAAGCGCGT'
    elif strain == 'bSHIVAD8EOM':
        inner_left =  'ACGCGCAGCGTA'
        inner_right = 'CAGCTAGCGCGT'
    elif strain == 'cSHIV174M':
        inner_left =  'ACCGGACAGCGT'
        inner_right = 'AGCTGACCCGGT'
    elif strain == 'cSHIV224M':
        inner_left =  'ACCGGACAGCGT'
        inner_right = 'CACAAGTCCGGT'
    else:
        raise ValueError(f"Unknown strain: {strain}")

    return {'barcode_primer': f'({left}{inner_left}){{i<={i},d<={d},s<={s}}}(?<umi>.{{10}})({inner_right}{right}){{i<={i},d<={d},s<={s}}}'}

# Built-in regex patterns
BUILTIN_PATTERNS = {
    'SIVNFL': {
        'left_primer': '(CAAGCAGAAGACGGCATACGAGAT){i<=2,d<=2,s<=4}(?<umi>.{18})(atatacttagaaaaggaagaaggcatcataccaga){i<=2,d<=2,s<=4}',
        'right_primer': '(gaagaactccgaaaaaggctaaggc){i<=2,d<=2,s<=4}(?<umi>.{18})(GATCTCGGTGGTCGCCGTATCATT){i<=2,d<=2,s<=4}',
    }
}


def get_regex_patterns(params):
    """Get regex patterns from either builtin or custom parameters"""
    
    if 'builtin' in params:
        patterns = {}
        for pattern in params['builtin'].split(','):
            if pattern not in BUILTIN_PATTERNS:
                patterns.update(_build_siv_barcode(pattern, 2, 2, 4))
            else:
                patterns.update(BUILTIN_PATTERNS[pattern])
    else:
        required = ['left_primer', 'right_primer', 'barcode_primer']
        missing = [p for p in required if p not in params]
        if missing:
            raise ValueError(f"Missing required primer patterns: {missing}")
        patterns = {
            'left_primer': params['left_primer'],
            'right_primer': params['right_primer'],
            'barcode_primer': params['barcode_primer']
        }
    
    return {
        'LEFT_PRIMER': regex.compile(patterns['left_primer'], flags=regex.BESTMATCH+regex.IGNORECASE),
        'RIGHT_PRIMER': regex.compile(patterns['right_primer'], flags=regex.BESTMATCH+regex.IGNORECASE),
        'BARCODE_PRIMER': regex.compile(patterns['barcode_primer'], flags=regex.BESTMATCH+regex.IGNORECASE)
    }

def extract_umi(reg, text, default=None):
    """Extract UMI sequence from text using regex pattern"""
    matches = reg.finditer(text)
    for match in matches:
        return match.capturesdict()['umi'][0]
    return default

@contextmanager
def make_new_bamfile(insam, outpath):
    """Context manager for creating new BAM file"""
    with open(outpath, mode='wb') as handle:
        header = insam.header.to_dict()
        
        # Add program group to header
        if 'PG' not in header:
            header['PG'] = []
            
        header['PG'].append({
            'ID': 'barcode_extract',
            'PN': 'damlab-barcode-extract',
            'VN': __version__,
            'DS': 'Extract barcodes and UMIs from sequence reads'
        })
        
        outsam = pysam.AlignmentFile(handle, mode='wb', header=header)
        yield outsam

def generate_metrics(reads_processed: dict) -> dict:
    """Generate metrics from barcode extraction results.
    
    Args:
        reads_processed: Dictionary tracking successful extractions
        
    Returns:
        Dictionary of metrics including total reads and success rates
    """
    metrics = {
        'total_reads': reads_processed['total'],
        'extraction_counts': {
            'barcode': reads_processed['barcode'],
            'left_umi': reads_processed['left_umi'],
            'right_umi': reads_processed['right_umi']
        },
        'pairwise_counts': {
            'barcode_left': reads_processed['barcode_left'],
            'barcode_right': reads_processed['barcode_right'],
            'left_right': reads_processed['left_right'],
            'all_three': reads_processed['all']
        }
    }
    
    return metrics

# Get input/output paths
bam_in_path = str(snakemake.input[0])
bam_out_path = str(snakemake.output[0])

# Get parameters with defaults
bc_tag = snakemake.params.get('barcode_tag', 'CR')  # Cellular barcode sequence
umi_tag = snakemake.params.get('umi_tag', 'OX')  # Original UMI bases
mapped_only = snakemake.params.get('mapped_only', False)  # Only process mapped reads
metrics_file = snakemake.output.get('metrics', None)  # Optional metrics output file
sample_name = snakemake.params.get('sample_name', None)  # Optional sample name

# Get regex patterns
patterns = get_regex_patterns(dict(snakemake.params))

# Track extraction statistics
reads_processed = {
    'total': 0,
    'mapped': 0,
    'barcode': 0,
    'left_umi': 0,
    'right_umi': 0,
    'barcode_left': 0,
    'barcode_right': 0,
    'left_right': 0,
    'all': 0
}

# Process BAM file
with open(bam_in_path, mode='rb') as handle:
    samfile = pysam.AlignmentFile(handle)
    with make_new_bamfile(samfile, bam_out_path) as newsam:
        for read in samfile:
            reads_processed['total'] += 1
            
            # Skip unmapped reads if mapped_only is True
            if mapped_only and read.is_unmapped:
                continue
            
            reads_processed['mapped'] += not read.is_unmapped

            if read.query_sequence:
                # Track successful extractions
                has_barcode = False
                has_left = False
                has_right = False
                
                # Extract barcode
                barcode = extract_umi(patterns['BARCODE_PRIMER'], read.query_sequence, default='N'*10)
                if not all(b == 'N' for b in barcode):
                    reads_processed['barcode'] += 1
                    has_barcode = True
                read.set_tag(bc_tag, barcode, 'Z')
                
                # Extract UMIs
                left_umi = extract_umi(patterns['LEFT_PRIMER'], read.query_sequence, default='N'*18)
                if not all(b == 'N' for b in left_umi):
                    reads_processed['left_umi'] += 1
                    has_left = True
                    
                right_umi = extract_umi(patterns['RIGHT_PRIMER'], read.query_sequence, default='N'*18)
                if not all(b == 'N' for b in right_umi):
                    reads_processed['right_umi'] += 1
                    has_right = True
                    
                read.set_tag(umi_tag, left_umi+right_umi, 'Z')
                
                # Track pairwise successes
                if has_barcode and has_left:
                    reads_processed['barcode_left'] += 1
                if has_barcode and has_right:
                    reads_processed['barcode_right'] += 1
                if has_left and has_right:
                    reads_processed['left_right'] += 1
                if has_barcode and has_left and has_right:
                    reads_processed['all'] += 1
                
                newsam.write(read)

# Write metrics if requested
if metrics_file:
    metrics = generate_metrics(reads_processed)
    if sample_name:
        metrics['sample_name'] = sample_name
    with open(metrics_file, 'w') as f:
        f.write('# Barcode extraction metrics\n')
        yaml.dump(metrics, f, default_flow_style=False)
