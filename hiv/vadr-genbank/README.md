# VADR GenBank Normaliser (`hiv/vadr-genbank`)

A Snakemake wrapper that strips NCBI version suffixes (`.1`, `.2`, …) from
GenBank flat files and FASTA files before they are passed to
[`hiv/vadr-build`](../vadr-build/).

## Why is this needed?

VADR's GenBank parser (`sqf_GenbankParse` in `sequip/sqp_seqfile.pm`) uses
the **first token of the LOCUS line** as its primary accession key, then
validates that the VERSION line (after stripping its own `.N` suffix) equals
that key exactly.

Two problems arise with files downloaded from NCBI:

**Problem 1 — legacy LOCUS name differs from the accession.**
`K03455.gb` has a historical LOCUS name `HIVHXB2CG` that pre-dates the modern
convention of LOCUS == accession:

```
LOCUS       HIVHXB2CG   9719 bp …   ← VADR uses this as $acc
VERSION     K03455.1                ← stripped to K03455 ≠ HIVHXB2CG → ERROR
```

**Problem 2 — VERSION contains a `.N` version suffix.**
Even when LOCUS == accession, VADR only strips the version suffix with its
own `seq_StripVersion` call *after* the LOCUS check.  If the LOCUS name is
already the bare accession this is fine, but when Problem 1 is also present
the `.N` suffix is a secondary concern.

Without normalisation `v-build.pl` aborts with:

```
ERROR in sqf_GenbankParse … version/accession mismatch for K03455,
line: VERSION K03455.1
```

This wrapper rewrites the affected fields in a new copy of the GenBank file;
the original file is never modified.

## Usage

```python
rule vadr_genbank:
    input:
        genbank = "refs/K03455.gb",
        fasta   = "refs/K03455.fasta",
    output:
        genbank = "refs/normalised/K03455.gb",
        fasta   = "refs/normalised/K03455.fasta",
    params:
        accession = "K03455",
    log:
        "logs/vadr_genbank_K03455.log"
    wrapper:
        "file:path/to/damlab-wrappers/hiv/vadr-genbank"

rule vadr_build:
    input:
        genbank = "refs/normalised/K03455.gb",
        fasta   = "refs/normalised/K03455.fasta",
    output:
        directory("vadr_models/K03455"),
    params:
        accession = "K03455",
        group     = "HIV",
        subgroup  = "HIV-1",
    log:
        "logs/vadr_build_K03455.log"
    wrapper:
        "file:path/to/damlab-wrappers/hiv/vadr-build"
```

## Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `accession` | **yes** | Bare accession (no `.N` suffix), e.g. `"K03455"` |
| `version` | no | Assert a specific wrapper version (warns if mismatch) |

## Inputs

| Name | Required | Description |
|------|----------|-------------|
| `genbank` | at least one | NCBI GenBank flat file (`.gb`) to normalise |
| `fasta` | at least one | NCBI FASTA file to normalise |

Either or both may be supplied.

## Outputs

| Name | Matches input | Description |
|------|--------------|-------------|
| `genbank` | `input.genbank` | Normalised GenBank file |
| `fasta` | `input.fasta` | Normalised FASTA file |

## What gets changed

**GenBank file** — two transformations are applied:

1. The LOCUS name token is replaced with the bare accession if it differs
   (e.g. `HIVHXB2CG` → `K03455`).
2. The `.N` version suffix is stripped from the VERSION line
   (e.g. `K03455.1` → `K03455`).

All other content is preserved verbatim.  If neither change is needed the
file is copied unchanged.

**FASTA file** — copied to `output.fasta` without modification.

> **Why is the FASTA not changed?**  `v-build.pl` requires FASTA headers to
> be in `<accession>.<version>` format (e.g. `>K03455.1`).  The Perl check
> unconditionally fails if the header has no `.N` suffix, so the version
> suffix **must** be preserved.  Only the GenBank VERSION line comparison uses
> the bare accession, which is why normalisation is limited to the `.gb` file.

## Custom references (no version suffix)

For custom reference sequences (e.g. a modified HIV genome with an inserted
reporter gene), you can use a plain identifier like `MY_HIV_GFP` that has no
version suffix.  In this case this wrapper is a no-op — the files are simply
copied — but it is still safe to include in the pipeline for consistency.

## Author

* Will Dampier, PhD

## Software Requirements

* Python ≥ 3.9 (standard library only — no additional packages required)
