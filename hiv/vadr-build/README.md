# VADR Build Wrapper (`hiv/vadr-build`)

A Snakemake wrapper for [`v-build.pl`](https://github.com/ncbi/vadr/blob/master/documentation/build.md),
which creates VADR homology model files from a GenBank/RefSeq accession.  The
resulting model directory can be passed to the `hiv/vadr-annotate` wrapper to
validate and annotate similar sequences.

> **Note:** `v-build.pl` fetches the accession from NCBI at runtime and
> requires an internet connection.  Only sequences ≤ 25 kb are supported.

## Usage

```python
rule vadr_build:
    output:
        directory("vadr_models/K03455")
    params:
        accession = "K03455",   # HIV-1 HXB2 RefSeq accession
        group    = "HIV",       # optional --group label
        subgroup = "HIV-1",     # optional --subgroup label
        extra    = "",          # any additional v-build.pl flags
    log:
        "logs/vadr_build_K03455.log"
    wrapper:
        "file:path/to/damlab-wrappers/hiv/vadr-build"
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `accession` | **yes** | — | GenBank / RefSeq accession to build a model from |
| `group` | no | `""` | Value passed to `--group` (e.g. `"HIV"`) |
| `subgroup` | no | `""` | Value passed to `--subgroup` (e.g. `"HIV-1"`) |
| `extra` | no | `""` | Arbitrary additional flags forwarded verbatim to `v-build.pl` |
| `version` | no | — | Assert a specific wrapper version (warns if mismatch) |

## Input

By default there are no file inputs — `v-build.pl` fetches the sequence and
feature table directly from NCBI.  Two optional named inputs enable **offline
operation**:

| Name | Description |
|------|-------------|
| `fasta` | Local FASTA file for the accession (passed as `--infa`). Skips the NCBI FASTA fetch. |
| `genbank` | Local GenBank flat file (passed as `--gb --ingb`). Skips the NCBI feature-table fetch and provides feature annotations from the local file. Typically used together with `fasta`. |

**Offline example:**

```python
rule vadr_build_offline:
    input:
        fasta   = "refs/K03455.fa",      # pre-normalised by hiv/vadr-genbank
        genbank = "refs/K03455.gb",      # pre-normalised by hiv/vadr-genbank
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

> **Note on local files from NCBI:** Files downloaded directly from NCBI
> (e.g. `K03455.gb`, `K03455.fasta`) contain versioned accessions such as
> `VERSION K03455.1` and `>K03455.1`.  VADR's GenBank parser compares these
> against the bare accession argument and fails with a *version/accession
> mismatch* error.  Use the `hiv/vadr-genbank` wrapper as a pre-processing
> step to strip the version suffixes before passing the files here.

## Output

A directory containing the VADR model files for the accession:

| File | Description |
|------|-------------|
| `<acc>.cm` | Infernal covariance model |
| `<acc>.minfo` | Model feature information |
| `<acc>.*.[np]hr/nin/nsq` | BLAST databases |
| `<acc>.fa` | Reference FASTA |
| `<acc>.stk` | Stockholm alignment |

Pass this directory as `input.mdir` to `hiv/vadr-annotate`.

## Error Handling

The wrapper raises errors for:
- Missing or empty `params.accession`
- Output directory already exists and is non-empty
- `v-build.pl` execution failures (propagated through Snakemake log)

## Author

* Will Dampier, PhD

## Software Requirements

* [VADR](https://github.com/ncbi/vadr) ≥ 1.6.4 (via `bioconda::vadr`)
* Internet access (for NCBI accession fetch)
