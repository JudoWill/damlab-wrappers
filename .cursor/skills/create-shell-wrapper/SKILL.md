---
name: damlab-shell-wrapper
description: >-
  Create Snakemake shell-command wrappers in damlab-wrappers. Use when adding
  a new wrapper for a CLI tool, scaffolding test files, wiring a wrapper into
  the root Snakefile, or following wrapper conventions in this repository.
  Triggers on: "new wrapper", "add wrapper", "create wrapper", "test snakefile",
  "wrapper boilerplate", "wrapper convention".
---

# damlab-wrappers: Shell Wrapper Conventions

## Directory structure

Every wrapper lives at `{category}/{tool-name}/`. Always read an existing
complete wrapper (e.g. `CRISPR/crispresso-core/`) before writing a new one.

```
{category}/{tool}/
├── wrapper.py
├── environment.yaml
├── README.md
└── test/
    ├── Snakefile
    ├── tests.py
    ├── env.yaml
    └── test_data/        # small synthetic inputs, committed to git
```

---

## `wrapper.py` template

```python
"""One-line description of what this wrapper does."""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

from snakemake.shell import shell  # type: ignore

if "snakemake" not in locals():
    import snakemake  # type: ignore

if hasattr(snakemake, "params") and "version" in snakemake.params:
    requested_version = snakemake.params.version
    if requested_version != __version__:
        print(f"Warning: Requested version {requested_version} does not match "
              f"wrapper version {__version__}")

# --- Inputs ---
input_file = snakemake.input[0]                         # positional
opt_input  = snakemake.input.get("named_key", None)     # optional named

# --- Outputs ---
output_file = snakemake.output[0]

# --- Parameters ---
required_param = snakemake.params.required_param         # raise AttributeError if missing
opt_param      = snakemake.params.get("opt_param", None)
flag_param     = snakemake.params.get("flag_param", False)
extra          = snakemake.params.get("extra", "")       # always include
threads        = snakemake.threads

# --- Build command ---
args = [
    f"--required {required_param}",
    f"--threads {threads}",
]
if opt_param is not None:
    args.append(f"--opt {opt_param}")
if opt_input:
    args.append(f"--named {opt_input}")
if flag_param:
    args.append("--flag")
if extra:
    args.append(extra)

args_str = " ".join(args)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(f"tool-command {args_str} {input_file} > {output_file} {log}")
```

**Parameter access rules:**
- Use `snakemake.params.get("key", default)` for optional params
- Use `snakemake.params.key` (no default) for required params — raises `AttributeError` with a clear name
- Use `snakemake.input.get("key", None)` for optional named inputs
- Always include `extra` as a verbatim pass-through escape hatch

---

## `environment.yaml` template

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - tool-name>=X.Y
  - snakemake-wrapper-utils>=0.3
```

---

## Test Snakefile conventions

**Naming:**
- Root variable: `TEST_ROOT_{CATEGORY}__{tool} = Path(workflow.snakefile).parent`
  - Double underscore separates category from tool; use `__` between all segments
  - e.g. `TEST_ROOT_CRISPR__crispresso_core`
- All-rule: `test_{CATEGORY}__{tool}__all`
- Individual rules: `{CATEGORY}__{tool}__test_{variant}`
- Wrapper reference: `f"file:{TEST_ROOT_{CATEGORY}__{tool}.parent}"`

**Structure of the `__all` rule:**
```python
rule test_{CATEGORY}__{tool}__all:
    input:
        files = [TEST_ROOT / "test_output/result_file"],   # plain paths, NOT directory()
        tests = [TEST_ROOT / "tests.py"],
    conda: "env.yaml"
    log:
        stdout = TEST_ROOT / "test_output/pytest.stdout.log",
        stderr = TEST_ROOT / "test_output/pytest.stderr.log",
    shell: "cd {TEST_ROOT} && pytest {input.tests} > {log.stdout} 2> {log.stderr}"
```

**Individual test rules:**
```python
rule {CATEGORY}__{tool}__test_{variant}:
    input:
        TEST_ROOT / "test_data/input_file",
    output:
        TEST_ROOT / "test_output/result_file",          # file output
        # OR:
        directory(TEST_ROOT / "test_output/result_dir"), # directory output
    params:
        param = value,
    threads: 1
    log:
        temp(TEST_ROOT / "test_output/variant.log"),    # temp() cleans up job logs
    wrapper:
        f"file:{TEST_ROOT.parent}"
```

**Test `env.yaml`** (minimal — just pytest):
```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - pytest
```

---

## Root Snakefile wiring

Add to [`Snakefile`](../../Snakefile) (two places):

```python
# Near other includes:
include: "{category}/{tool}/test/Snakefile"

# Inside rule test_cpu input list:
rules.test_{CATEGORY}__{tool}__all.log,
```

---

## Critical gotchas

### `directory()` only on outputs
Never use `directory()` in `input:` blocks — Snakemake raises an error.
When referencing a `directory()` output from another rule, use a plain path.

### `cd` + relative paths
If the wrapper must `cd` before calling the tool (e.g. tools that write to cwd),
all paths must be made absolute **before** the `cd`:

```python
import os
from pathlib import Path

parent_dir = str(Path(snakemake.output[0]).resolve().parent)
input_dirs = [os.path.abspath(d) for d in input_dirs]

# log_fmt_shell() returns relative paths — they break after cd
# Build the log redirect manually:
if snakemake.log:
    log = f"> {os.path.abspath(str(snakemake.log[0]))} 2>&1"
else:
    log = ""

shell(f"cd {parent_dir} && tool {args_str} {log}")
```

### Tools that write `Tool_on_{name}/` output directories
Derive `--name` from the output directory basename and pass the parent as `--output_folder`:

```python
out_path = Path(snakemake.output[0])
bname = out_path.name
name = snakemake.params.get("name", None)
if name is None:
    prefix = "Tool_on_"
    if bname.startswith(prefix):
        name = bname[len(prefix):]
    else:
        raise ValueError(
            f"Output directory '{bname}' must start with '{prefix}' "
            "or params.name must be set explicitly."
        )
args.append(f"--output_folder {out_path.parent}")
args.append(f"--name {name}")
```

### Variable-length input lists
For wrappers that accept any number of inputs:

```python
if hasattr(snakemake.input, "named_key"):
    items = snakemake.input.named_key
    if isinstance(items, str):   # single input comes back as str
        items = [items]
else:
    items = list(snakemake.input)
```

### Non-repetitive test sequences
For bioinformatics tools with windowed alignment (e.g. CRISPResso): if a guide
or probe sequence appears at many positions in a repetitive amplicon, tools may
fail with quantification window errors. Always use unique, non-repetitive
sequences in test data.
