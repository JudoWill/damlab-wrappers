"""Wrapper for tree rerooting using DendroPy.

This wrapper provides a Snakemake interface to DendroPy's tree rerooting functionality,
allowing users to reroot phylogenetic trees at a specified taxon. This is particularly
useful for ensuring proper tree orientation and outgroup placement.
"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.1"

import os
from typing import Optional

import dendropy  # type: ignore

# This is a common pattern in Snakemake wrappers
# It allows the wrapper to be imported without snakemake being in the global namespace
# This is useful for testing and linting
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Check if version is specified and compatible
if hasattr(snakemake, "params") and "version" in snakemake.params:
    requested_version = snakemake.params.version
    if requested_version != __version__:
        print(f"Warning: Requested version {requested_version} does not match wrapper version {__version__}")

# Get input and output files
input_tree: str = snakemake.input[0]
output_tree: str = snakemake.output[0]

# Validate input file exists
if not os.path.exists(input_tree):
    raise FileNotFoundError(f"Input file {input_tree} does not exist")

# Get parameters with defaults
# This is a common pattern in wrappers - providing sensible defaults
root_taxon: Optional[str] = snakemake.params.get("root_taxon")
schema: str = snakemake.params.get("schema", "newick")  # Default to newick format
preserve_branch_lengths: bool = snakemake.params.get("preserve_branch_lengths", True)
preserve_support_values: bool = snakemake.params.get("preserve_support_values", True)

# Validate required parameters
if not root_taxon:
    raise ValueError("root_taxon parameter must be specified")

# Read the tree
try:
    tree = dendropy.Tree.get(
        path=input_tree,
        schema=schema
    )
except Exception as e:
    raise ValueError(f"Failed to read tree file: {str(e)}")

# Find the node corresponding to the root taxon
root_node = None
for node in tree.leaf_node_iter():
    if node.taxon and node.taxon.label == root_taxon:
        root_node = node
        break

if not root_node:
    raise ValueError(f"Root taxon '{root_taxon}' not found in tree")

# Reroot the tree
try:
    tree.reroot_at_node(
        root_node,
        update_bipartitions=True
    )
except Exception as e:
    raise ValueError(f"Failed to reroot tree: {str(e)}")

# Write the rerooted tree
try:
    tree.write(
        path=output_tree,
        schema=schema,
        suppress_edge_lengths=not preserve_branch_lengths,
        suppress_internal_node_labels=not preserve_support_values
    )
except Exception as e:
    raise ValueError(f"Failed to write rerooted tree: {str(e)}") 