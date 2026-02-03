import os
from pathlib import Path
from snakemake.utils import min_version

# Require Snakemake 7.0 or higher for better module support
min_version("7.0")

# Get the base directory of the project
BASE_DIR = Path(workflow.basedir)
WRAPPER_VERSION = 'v6.0.0'

# Common configuration
# configfile: "test_configs/default.txt"

# Include all test Snakefiles
# These will be loaded dynamically based on the test configuration
def get_test_snakefiles():
    with open(config["test_config"], "r") as f:
        tools = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    return [f"{tool}/test/Snakefile" for tool in tools]

# TODO: Once migrated, this will do all tests
# Include all test Snakefiles
# for snakefile in get_test_snakefiles():
#    include: snakefile

include: "example/shell/test/Snakefile"
include: "example/api/test/Snakefile"
include: "barcode/extract/test/Snakefile"
include: "barcode/correct/test/Snakefile"
include: "cigarmath/bam2gel/test/Snakefile"
include: "cigarmath/deletion_frequency/test/Snakefile"
include: "cigarmath/maximal_deletion_size/test/Snakefile"
include: "cigarmath/pileup/test/Snakefile"
include: "cigarmath/slice/test/Snakefile"

include: "dorado/demux/test/Snakefile"
include: "dorado/duplex/test/Snakefile"
include: "dorado/simplex/test/Snakefile"

include: "huggingface/hiv-bert/test/Snakefile"

include: "MSA/muscle/test/Snakefile"
include: "pandas/merge/test/Snakefile"

include: "phylo/FastTree/test/Snakefile"
include: "phylo/phytreeviz/test/Snakefile"
include: "phylo/reroot/test/Snakefile"

# Rule to run all tests with coverage


rule test_cpu:
    input:
        # From each test, capture the final log file   
        # With --forceall this will rerun the entire suite
        rules.test_example__shell__all.log,
        rules.test_example__api__all.log,
        rules.test_barcode__extract__all.log,
        rules.test_barcode__correct__all.log,
        rules.test_cigarmath__bam2gel__all.log,
        rules.test_cigarmath__DF__all.log,
        rules.test_cigarmath__MDS__all.log,
        rules.test_cigarmath__pileup__all.log,
        rules.test_cigarmath__slice__all.log,
        rules.test_dorado__demux__all.log,
        rules.test_msa__muscle__all.log,
        rules.test_pandas__merge__all.log,
        rules.test_phylo__fasttree__all.log,
        rules.test_phylo__phytreeviz__all.log,
        rules.test_phylo__reroot__all.log


rule test_gpu:
    input:
        rules.test_dorado__duplex__all.log,
        rules.test_dorado__simplex__all.log,
        rules.test_huggingface__hivbert__all.log


rule test_all:
    input:
        rules.test_cpu.input,
        rules.test_gpu.input,


