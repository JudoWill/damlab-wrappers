# VADR Annotate Wrapper (`hiv/vadr-annotate`)

A Snakemake wrapper for [`v-annotate.pl`](https://github.com/ncbi/vadr/blob/master/documentation/annotate.md),
which classifies and annotates viral sequences using VADR homology models.
Sequences are reported as *passing* (zero fatal alerts) or *failing* (one or
more fatal alerts), and feature tables / alignment details are written to the
output directory.

Build your model library first with the `hiv/vadr-build` wrapper, or point
`input.mdir` at an existing VADR model library.

## Usage

```python
rule vadr_annotate:
    input:
        sequences = "clipped/{sample}.fa",
        mdir      = "vadr_models/K03455",   # optional — omit to use VADR defaults
    output:
        directory("vadr/{sample}"),
    params:
        mode         = "hiv",  # preset --alt_pass for all recoverable HIV alerts
        noseqnamemax = True,   # suppress long-name warnings (common for haplotype IDs)
        extra        = "",     # any additional v-annotate.pl flags
    threads: 4
    log:
        "logs/{sample}.vadr_annotate.log"
    wrapper:
        "file:path/to/damlab-wrappers/hiv/vadr-annotate"
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `mode` | no | `""` | Preset mode; `"hiv"` enables `--alt_pass` for all recoverable alerts (see below) |
| `noseqnamemax` | no | `False` | Pass `--noseqnamemax`; suppresses errors for long sequence names |
| `extra` | no | `""` | Arbitrary additional flags forwarded verbatim to `v-annotate.pl` |
| `version` | no | — | Assert a specific wrapper version (warns if mismatch) |

### `mode="hiv"` — HIV-aware alert thresholds

VADR defines many alerts as *fatal by default*, which causes HIV sequences with
typical proviral variation (frameshifts, internal stop codons, deletion, etc.)
to fail annotation even when the sequence is biologically valid.

Setting `mode="hiv"` passes `--alt_pass` with every alert code that is *fatal
by default but not always-fatal*, effectively converting them all to
non-fatal warnings.  The only alert that can never be overridden is `ftskipfl`
(always-fatal).

Alerts downgraded by `mode="hiv"`:

| Category | Alert codes |
|----------|-------------|
| Sequence-level | `incsbgrp`, `incgroup`, `lowcovrg`, `dupregin`, `discontn`, `indfstrn`, `lowsim5s`, `lowsim3s`, `lowsimis`, `nmiscftr`, `deletins` |
| Feature-level (CDS) | `mutstart`, `mutendcd`, `mutendns`, `mutendex`, `unexleng`, `cdsstopn`, `cdsstopp`, `fsthicft`, `fsthicfi`, `fstukcft`, `fstukcfi`, `mutspst5`, `mutspst3`, `peptrans`, `pepadjcy` |
| Feature-level (annotation) | `indfantp`, `indfantn`, `indf5gap`, `indf5lcn`, `indf5plg`, `indf5pst`, `indf3gap`, `indf3lcn`, `indf3plg`, `indf3pst`, `indfstrp`, `insertnp`, `deletinp`, `deletinf` |
| Feature-level (similarity) | `lowsim5n`, `lowsim5l`, `lowsim3n`, `lowsim3l`, `lowsimin`, `lowsimil` |

You can further customise alert behaviour via `params.extra`, e.g.:
```python
extra = "--alt_fail lowcovrg"  # re-promote specific alerts back to fatal
```

## Input

| Name | Required | Description |
|------|----------|-------------|
| `sequences` | **yes** | FASTA file of sequences to annotate |
| `mdir` | no | Path to a VADR model library directory (from `vadr-build`). If omitted, VADR's built-in default library is used |

## Output

A directory containing all `v-annotate.pl` output files:

| File pattern | Description |
|-------------|-------------|
| `*.pass.fa` | Sequences with zero fatal alerts |
| `*.fail.fa` | Sequences with ≥ 1 fatal alert |
| `*.sqc` | Per-sequence classification summary |
| `*.ftr` | Per-feature annotation table |
| `*.alt` | Alert detail listing |
| `*.log` | Full execution log |

## Threading

When `threads > 1`, `--split --cpu <threads>` is passed automatically to
enable multi-threaded processing.

## Error Handling

The wrapper raises errors for:
- Missing `input.sequences` file
- Non-existent `input.mdir` directory
- Output directory already exists and is non-empty
- `v-annotate.pl` execution failures (propagated through Snakemake log)

## Author

* Will Dampier, PhD

## Software Requirements

* [VADR](https://github.com/ncbi/vadr) ≥ 1.6.4 (via `bioconda::vadr`)
