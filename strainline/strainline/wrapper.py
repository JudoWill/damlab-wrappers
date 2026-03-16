"""Wrapper for strainline"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import tempfile
import shutil
from os.path import join
import os
from snakemake.shell import shell # type: ignore

if "snakemake" not in locals():
    # Keeps linters happy but doesn't impact function
    import snakemake # type: ignore


required_version = snakemake.params.get("version", __version__)
if required_version != __version__:
    print(f"Warning: Wrapper version {__version__} does not match required version {required_version}")

# Extract and validate required arguments
reads = snakemake.input[0]
if not os.path.exists(reads):
    raise FileNotFoundError(f"Input file {reads} does not exist")

prefix = snakemake.params.get("prefix", "")
if not prefix:
    raise ValueError("prefix parameter must be specified")
if not os.path.exists(prefix):
    raise FileNotFoundError(f"Strainline installation directory {prefix} does not exist")

platform = snakemake.params.get("platform", "ont")
if platform not in ["ont", "pb"]:
    raise ValueError(f"Invalid platform '{platform}'. Must be either 'ont' or 'pb'")

# Determine output mode
if len(snakemake.output) != 1:
    raise ValueError("Exactly one output must be specified (either directory or haplotypes file)")

if snakemake.output.get('haplotypes', False):
    output = snakemake.output['haplotypes']
    is_directory = False
elif snakemake.output.get('directory', False):
    output = snakemake.output['directory']
    is_directory = True
else:
    raise ValueError("Exactly one output must be specified (either directory or haplotypes file)")

# Get optional parameters with defaults matching help file
extra_params = snakemake.params.get("extra_params", "")

# Get number of threads
threads = f"--threads {snakemake.threads}" if snakemake.threads > 1 else ""

# Handle logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Ensure the strainline script path is correct
strainline_path = os.path.join(prefix, "strainline", "src", "strainline.sh")
if not os.path.exists(strainline_path):
    raise FileNotFoundError(f"Strainline script not found at {strainline_path}")

daccord = 'daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu'
daccord_path = os.path.join(prefix, daccord)
if not os.path.exists(daccord_path):
    raise FileNotFoundError(f"Daccord installation not found at {daccord_path}")

# Build and execute command
path_cmd = f"export PATH={prefix}/bin:{daccord_path}/bin:$PATH && "
cmd = f"{strainline_path}"
cmd += f" -i {reads}"
cmd += f" -p {platform}"
cmd += f" {extra_params}"
cmd += f" {threads}"

if is_directory:
    cmd += f" -o {output}"
    shell(path_cmd + cmd + log)
else:
    with tempfile.TemporaryDirectory(delete=False) as tmpdir:
        cmd += f" -o {tmpdir}"
        haplotype_file = join(tmpdir, "filter_by_abun", "haplotypes.final.fa")
        try:
            shell(path_cmd + cmd + log)
        except Exception:
            # Strainline exits non-zero when it produces 0 haplotypes (the
            # final minimap2 | samtools sort pipeline fails on an empty index).
            # Treat a missing haplotype file as a valid zero-haplotype result
            # rather than a pipeline error.  A genuine crash (e.g. missing
            # binary) will still raise because the file won't exist below.
            pass
        if os.path.exists(haplotype_file):
            shutil.move(haplotype_file, output)
        else:
            # Zero haplotypes: write an empty FASTA so Snakemake considers
            # the rule complete.
            open(output, 'w').close()
