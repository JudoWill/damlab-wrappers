# VADR tbl2gbk Wrapper (`hiv/vadr-tbl2gbk`)

A Snakemake wrapper for [NCBI `table2asn`](https://www.ncbi.nlm.nih.gov/genbank/table2asn/),
which converts the NCBI 5-column feature table and passing-sequence FASTA
produced by `v-annotate.pl` into an annotated GenBank flat file (`.gbf`) and
an ASN.1 submission file (`.sqn`).

This wrapper is designed to be used immediately after `hiv/vadr-annotate`,
consuming its `.vadr.pass.fa` and `.vadr.pass.tbl` outputs.

## VADR output → GenBank pathway

```
vadr/{sample}/
  {sample}.vadr.pass.fa    ─┐
  {sample}.vadr.pass.tbl   ─┴─→  table2asn  →  {sample}.gbf
                                               {sample}.sqn  (optional)
```

## Usage

### Preferred: directory mode (use after `hiv/vadr-annotate`)

Pass the `vadr-annotate` output directory directly.  The wrapper locates
`*.vadr.pass.fa` and `*.vadr.pass.tbl` inside it automatically.  This avoids
a `MissingInputException` that would occur if files inside a Snakemake
`directory()` output were referenced as explicit inputs.

```python
rule vadr_tbl2gbk:
    input:
        vadr_dir = "vadr/{sample_name}",   # directory() output of vadr-annotate
    output:
        gbf = "gbk/{sample_name}.gbf",
        sqn = "gbk/{sample_name}.sqn",     # optional — omit if not needed
    params:
        organism = "Human immunodeficiency virus 1",
        taxid    = 11676,
        extra    = "",
    log:
        "logs/{sample_name}.tbl2gbk.log"
    wrapper:
        "file:path/to/damlab-wrappers/hiv/vadr-tbl2gbk"
```

### Alternative: explicit file mode

```python
rule vadr_tbl2gbk:
    input:
        sequences     = "vadr/{sample_name}/{sample_name}.vadr.pass.fa",
        feature_table = "vadr/{sample_name}/{sample_name}.vadr.pass.tbl",
    output:
        gbf = "gbk/{sample_name}.gbf",
    params:
        organism = "Human immunodeficiency virus 1",
    log:
        "logs/{sample_name}.tbl2gbk.log"
    wrapper:
        "file:path/to/damlab-wrappers/hiv/vadr-tbl2gbk"
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `organism` | no | `""` | Organism name inserted as a source feature qualifier (`-j "[organism=...]"`) |
| `taxid` | no | `""` | NCBI taxonomy ID passed to `table2asn -T` for taxonomy lookup |
| `extra` | no | `""` | Arbitrary additional flags forwarded verbatim to `table2asn` |
| `version` | no | — | Assert a specific wrapper version (warns if mismatch) |

## Input

Two modes are supported — provide exactly one of the following:

| Mode | Input name | Description |
|------|-----------|-------------|
| **Directory** (preferred) | `vadr_dir` | `vadr-annotate` output directory; `*.vadr.pass.fa` and `*.vadr.pass.tbl` are located automatically |
| **Explicit** | `sequences` + `feature_table` | Pass the FASTA and `.tbl` files directly |

## Output

| Name | Required | Description |
|------|----------|-------------|
| `gbf` | **yes** (first output) | GenBank flat file (`.gbf`) — human-readable annotated sequences |
| `sqn` | no | ASN.1 submission file (`.sqn`) — required for NCBI submission |

If `output.sqn` is not declared, the `.sqn` file produced by `table2asn` is discarded.

## Notes

- The wrapper raises a warning (not an error) if the input feature table is
  empty; this occurs when all sequences failed VADR checks.  An empty `.gbf`
  is still written to satisfy the Snakemake DAG.
- `table2asn` is run in a temporary directory; inputs are symlinked under a
  fixed stem (`seqs`) to give predictable output file names before moving
  them to the requested output paths.
- Use `-V b` (already applied by default) to enable GenBank flat file
  generation alongside the `.sqn` file.

## Author

* Will Dampier, PhD

## Software Requirements

* [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/) ≥ 1.28 (via `bioconda::table2asn`)
