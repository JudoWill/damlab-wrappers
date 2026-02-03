"""Wrapper for Biopython translation"""

__author__ = "Example Author"
__copyright__ = "Copyright 2024"
__email__ = "example@example.com"
__license__ = "MIT"
__version__ = "1.1.0"

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import logging

# Configure logging to use Snakemake's log file
if "snakemake" in locals():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=snakemake.log[0] if hasattr(snakemake, 'log') and snakemake.log else None
    )
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

logger = logging.getLogger('example-api-wrapper')

# This is a common pattern in Snakemake wrappers
# It allows the wrapper to be imported without snakemake being in the global namespace
# This is useful for testing and linting
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Check if version is specified and compatible
if hasattr(snakemake, "params") and "version" in snakemake.params:
    requested_version = snakemake.params.version
    if requested_version != __version__:
        logging.warning(f"Requested version {requested_version} does not match wrapper version {__version__}")

# Get input and output files
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Log useful parameters
logger.info(f"Input file: {input_file}")
logger.info(f"Output file: {output_file}")

# Get parameters with defaults
# This is a common pattern in wrappers - providing sensible defaults
frame = snakemake.params.get("frame", 0)
table = snakemake.params.get("table", 1)
stop_symbol = snakemake.params.get("stop_symbol", "*")
to_stop = snakemake.params.get("to_stop", False)
cds = snakemake.params.get("cds", False)

# Get the translation table
# This is a common pattern in API wrappers - converting parameters to API objects
translation_table = CodonTable.ambiguous_dna_by_id[table]

# Process the sequences
# This is the core of an API wrapper - using the API to process the data
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    # Read input sequences
    for num, record in enumerate(SeqIO.parse(infile, "fasta")):
        # Get the sequence
        seq = record.seq
        
        # Apply frame if specified
        if frame > 0:
            seq = seq[frame:]
        
        # Translate the sequence
        # This is where we use the Biopython API
        translated = seq.translate(
            table=translation_table,
            stop_symbol=stop_symbol,
            to_stop=to_stop,
            cds=cds
        )
        
        # Create a new record with the translated sequence
        translated_record = record
        translated_record.seq = translated
        
        # Write the translated sequence
        SeqIO.write(translated_record, outfile, "fasta") 

logger.info(f"Processed {num + 1} sequences")