"""Run samtools consensus to generate a FASTA consensus sequence from a sorted BAM."""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026, Will Dampier"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

from snakemake.shell import shell  # type: ignore

if "snakemake" not in locals():
    import snakemake  # type: ignore

if hasattr(snakemake, "params") and "version" in snakemake.params:
    requested_version = snakemake.params.version
    if requested_version != __version__:
        print(
            f"Warning: Requested version {requested_version} does not match "
            f"wrapper version {__version__}"
        )

from snakemake_wrapper_utils.samtools import get_samtools_opts  # type: ignore

samtools_opts = get_samtools_opts(
    snakemake, parse_write_index=False, parse_output_format=False
)
mode = snakemake.params.get("mode", "bayesian")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "samtools consensus -f FASTA -m {mode} {samtools_opts} {extra}"
    " -o {snakemake.output[0]} {snakemake.input[0]} {log}"
)
