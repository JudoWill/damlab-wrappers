# BAM to CSV Wrapper

This wrapper converts BAM/SAM files to CSV format, allowing you to extract specific fields from each alignment.

## Input

- BAM/SAM file

## Output

- CSV file containing the specified fields for each alignment

## Parameters

- `fields` (list): List of field names to extract from each alignment. Can include:
  - Standard BAM fields (e.g., `query_name`, `query_sequence`, `reference_name`, `reference_start`)
  - BAM tags (e.g., `NM`, `MD`, `AS`)

## Example Usage

```python
rule bam2csv:
    input:
        "input.bam"
    output:
        "output.csv"
    params:
        fields=["query_name", "reference_name", "reference_start", "NM", "MD"]
    wrapper:
        "damlab-wrappers/cigarmath/bam2gel"
```

## Default Fields

If no fields are specified, the wrapper will extract:
- `query_name`
- `reference_name`
- `reference_start`

## Notes

- The wrapper handles both standard BAM fields and BAM tags
- If a field is not found in an alignment, it will be represented as `None` in the CSV
- The output CSV will have a header row containing the field names 