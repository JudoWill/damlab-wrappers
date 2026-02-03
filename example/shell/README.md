# SeqKit Translate Wrapper

Version: 1.0.0

This wrapper demonstrates how to create a Snakemake wrapper for a shell command. It wraps the `seqkit translate` command, which translates DNA/RNA sequences to protein sequences.

## Input
* FASTA/FASTQ file containing DNA/RNA sequences

## Output
* FASTA file containing translated protein sequences

## Parameters
* `frame` (optional, default: 1)
    Reading frame for translation (1, 2, 3, -1, -2, -3, or 6 for all six frames)
* `min_len` (optional, default: None)
    Minimum length of amino acid sequence
* `transl_table` (optional, default: 1)
    Translation table/genetic code number
* `trim` (optional, default: false)
    Whether to remove all 'X' and '*' characters from the right end of the translation
* `extra` (optional, default: "")
    Additional parameters to pass to seqkit translate

## Example
```python
rule translate_sequences:
    input:
        "sequences.fasta"
    output:
        "proteins.fasta"
    threads: 4
    params:
        frame=1,
        min_len=100,
        transl_table=1,
        trim=True,
        extra="--clean"
    wrapper:
        "file:path/to/damlab-wrappers/example/shell"
```

## Output Format
The output is a FASTA file containing the translated protein sequences.

## Author
* Example Author

## Software Requirements
* [seqkit](https://bioinf.shenwei.me/seqkit/) (tested with v0.16.1) 