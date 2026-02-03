# Regex Pattern Matching Wrapper

This wrapper searches for regex patterns in FASTA/FASTQ files and outputs the matches to a CSV file.

## Input

- FASTA/FASTQ file containing sequences to search (supports .fasta, .fa, .fna, .fastq, .fq, and their gzipped versions)
- Dictionary of regex patterns to search for

## Output

- CSV file containing:
  - Read name
  - One column per pattern containing the first match found (if any)
- Optional metrics file in YAML format

## Parameters

- `patterns`: Dictionary where:
  - Keys are the column names for the matches in the output CSV
  - Values are the regex patterns to search for
  - Each pattern will be compiled with BESTMATCH and IGNORECASE flags
  - The first match for each pattern will be recorded
  - **Note**: When using regex patterns with curly braces (e.g., `{5}`) in Snakemake, you must:
    1. Double the curly braces (`{{5}}`) to escape them from Snakemake's wildcard expansion
    2. Wrap the patterns in a lambda function to prevent early evaluation
    ```python
    params:
        patterns=lambda wildcards: {
            'start_pattern': r"ATG[ATGC]{{5}}TAG",  # Correct: double curly braces
            'end_pattern': r"GAT[ATGC]{{5}}CTA"     # Correct: double curly braces
        }
    ```
- `both_strands`: Boolean (default: False)
  - If True, searches both the forward and reverse complement strands
  - If a pattern is not found in the forward strand, the reverse complement is searched
  - The first match found (either forward or reverse) is returned

## Example Usage

```python
rule regex_match:
    input:
        seq="input.fasta"
    output:
        csv="matches.csv",
        metrics="metrics.yaml"
    params:
        patterns=lambda wildcards: {
            'start_pattern': r"ATG[ATGC]{{10,20}}TAG",  # Example pattern 1
            'end_pattern': r"GAT[ATGC]{{5,15}}CTA"      # Example pattern 2
        },
        both_strands=True  # Search both strands
    wrapper:
        "damlab-wrappers/regex/match"
```

## Output Format

The CSV output will have the following columns:
- `read_name`: The name of the read from the FASTA/FASTQ file
- One column for each pattern key in the patterns dictionary, containing the first match found (if any)

## Metrics

If a metrics file is specified, it will contain:
- Total number of reads processed
- Number of matches found for each pattern (using the pattern keys as names) 