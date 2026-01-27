from os.path import join
from snakemake.utils import validate

# Schema validation using Snakemake's built-in validate function
# Validate run metadata against nanopore_run schema

validate(config, join(WORKFLOW_DIR, "schemas", "nanopore_run.yaml"))


# Validate basecalling config against nanopore_basecalling schema
validate(config, join(WORKFLOW_DIR, "schemas", "nanopore_basecalling.yaml"))


# Validate samples DataFrame against nanopore_sample schema
validate(SAMPLES, join(WORKFLOW_DIR, "schemas", "nanopore_sample.yaml"))
