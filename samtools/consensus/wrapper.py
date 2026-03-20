"""Run samtools consensus to generate a FASTA consensus sequence from a sorted BAM."""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026, Will Dampier"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.1.0"

import os
import shlex
import subprocess
import tempfile
from pathlib import Path

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


def _parse_max_reads(value):
    """Return positive int cap, or None if downsampling is disabled."""
    if value is None or value is False:
        return None
    try:
        v = int(value)
    except (TypeError, ValueError):
        return None
    return v if v > 0 else None


def _job_tmpdir():
    r = getattr(snakemake, "resources", None)
    t = None
    if r is not None:
        if isinstance(r, dict):
            t = r.get("tmpdir")
        else:
            t = getattr(r, "tmpdir", None)
    return t or os.environ.get("TMPDIR") or "/tmp"


def _primary_log_path():
    lg = snakemake.log
    if isinstance(lg, str):
        return lg
    return lg[0]


def _count_alignments(bam_path: str, filter_flags: str, threads: int, log_fp) -> int:
    cmd = ["samtools", "view"]
    if threads and threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.extend(["-c", "-F", str(filter_flags), bam_path])
    proc = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    if proc.stderr:
        log_fp.write(proc.stderr)
    return int(proc.stdout.strip())


def _subsample_bam(
    bam_path: str,
    out_bam: str,
    fraction: float,
    seed: int,
    threads: int,
    log_fp,
) -> None:
    cmd = ["samtools", "view"]
    if threads and threads > 1:
        cmd.extend(["-@", str(threads)])
    # HTSlib: --subsample FLOAT, --subsample-seed INT (not --subsampling-seed; -s is INT.FRAC shorthand)
    cmd.extend(
        [
            "-b",
            "--subsample-seed",
            str(seed),
            "--subsample",
            str(fraction),
            "-o",
            out_bam,
            bam_path,
        ]
    )
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.stderr:
        log_fp.write(proc.stderr)
    proc.check_returncode()


samtools_opts = get_samtools_opts(
    snakemake, parse_write_index=False, parse_output_format=False
)
mode = snakemake.params.get("mode", "bayesian")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

max_reads = _parse_max_reads(snakemake.params.get("max_reads"))
bam_in = str(snakemake.input[0])
out_fa = str(snakemake.output[0])
threads = int(snakemake.threads) if snakemake.threads else 1

if max_reads is None:
    shell(
        "samtools consensus -f FASTA -m {mode} {samtools_opts} {extra}"
        " -o {snakemake.output[0]} {snakemake.input[0]} {log}"
    )
else:
    filter_flags = snakemake.params.get("count_filter_flags", "0x900")
    seed = int(snakemake.params.get("subsample_seed", 0))
    log_path = _primary_log_path()
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)

    with open(log_path, "a", encoding="utf-8") as logf:
        n = _count_alignments(bam_in, filter_flags, threads, logf)
        logf.write(
            f"[samtools/consensus wrapper] counted alignments (-F {filter_flags}): {n}; "
            f"max_reads={max_reads}\n"
        )
        logf.flush()

        consensus_input = bam_in
        tmp_bam = None
        try:
            if n > max_reads:
                frac = min(1.0, max_reads / n)
                logf.write(
                    f"[samtools/consensus wrapper] subsampling with --subsample {frac} "
                    f"--subsample-seed {seed}\n"
                )
                logf.flush()
                scratch = _job_tmpdir()
                Path(scratch).mkdir(parents=True, exist_ok=True)
                fd, tmp_bam = tempfile.mkstemp(suffix=".bam", dir=scratch)
                os.close(fd)
                _subsample_bam(bam_in, tmp_bam, frac, seed, threads, logf)
                consensus_input = tmp_bam

            q_in = shlex.quote(consensus_input)
            q_out = shlex.quote(out_fa)
            shell(
                "samtools consensus -f FASTA -m {mode} {samtools_opts} {extra}"
                f" -o {q_out} {q_in} {log}"
            )
        finally:
            if tmp_bam is not None and os.path.isfile(tmp_bam):
                os.unlink(tmp_bam)
