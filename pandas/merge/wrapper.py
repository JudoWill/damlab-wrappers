"""Wrapper for pandas merge functionality"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.1.2"

import pandas as pd
import yaml
import logging
from typing import List, Optional, Union, Tuple
from pathlib import Path

# Configure logging to use Snakemake's log file
if "snakemake" in locals():
    filename = snakemake.log[0] if hasattr(snakemake, 'log') and snakemake.log else None
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=filename
    )
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

logger = logging.getLogger('pandas-merge-wrapper')

if "snakemake" not in locals():
    import snakemake  # type: ignore

def merge_csv_files(input_files: List[str], 
                   on: Optional[Union[str, List[str], List[List[str]]]] = None,
                   how: str = 'inner',
                   suffixes: Optional[tuple] = None) -> pd.DataFrame:
    """Merge multiple CSV files iteratively.
    
    Args:
        input_files: List of CSV file paths to merge
        on: Column(s) to merge on. Can be:
            - A single column name (str) - used for all merges
            - A list of column names for all tables (List[str])
            - A list of lists where each inner list contains a key for that table (List[List[str]])
              In this case, the keys are used in a pairwise progressive pattern:
              For tables A, B, C with keys [['a'], ['b'], ['c']], the merges are:
              A.merge(B, left_on=['a'], right_on=['b'])
              result.merge(C, left_on=['b'], right_on=['c'])
        how: Type of merge to perform ('left', 'right', 'outer', 'inner')
        suffixes: Suffixes to apply to overlapping column names
        
    Returns:
        Merged DataFrame
    """
    if not input_files:
        raise ValueError("No input files provided")
        
    logger.info(f"Merging {len(input_files)} CSV files")
    logger.info(f"Merge parameters - on: {on}, how: {how}")
    
    # Read first file as base
    result = pd.read_csv(input_files[0])
    logger.info(f"Read base file: {input_files[0]} with shape {result.shape}")
    
    # If suffixes not provided, create default suffixes based on number of files
    if suffixes is None:
        suffixes = tuple(f"_{i}" for i in range(len(input_files)))
    
    # Enforce that the number of suffixes matches the number of files
    assert len(suffixes) == len(input_files), f"Number of suffixes ({len(suffixes)}) must match number of input files ({len(input_files)})"
    
    # Handle different types of 'on' parameter
    if isinstance(on, list) and on and isinstance(on[0], list):
        # Check if we have the right number of merge keys
        if len(on) != len(input_files):
            # If we have more merge keys than needed, use only the first ones
            if len(on) > len(input_files):
                logger.warning(f"Number of merge key lists ({len(on)}) is greater than number of input files ({len(input_files)}). Using only the first {len(input_files)} lists.")
                merge_keys = on[:len(input_files)]
            else:
                # If we have fewer merge keys than needed, repeat the last one
                logger.warning(f"Number of merge key lists ({len(on)}) is less than number of input files ({len(input_files)}). Repeating the last merge key list.")
                merge_keys = on + [on[-1]] * (len(input_files) - len(on))
        else:
            merge_keys = on
    else:
        # If on is a string or list of strings, use the same for all merges
        merge_keys = [on] * len(input_files)
    
    # Track the original key column names
    original_keys = []
    if isinstance(on, str):
        original_keys = [on] * len(input_files)
    elif isinstance(on, list) and on and isinstance(on[0], list):
        original_keys = [keys[0] if isinstance(keys, list) and len(keys) == 1 else None for keys in on]
    
    # Iteratively merge remaining files
    for i, file in enumerate(input_files[1:], 1):
        df = pd.read_csv(file)
        logger.info(f"Merging file: {file} with shape {df.shape}")
        
        # Get merge parameters for this iteration
        left_key = merge_keys[i-1]
        right_key = merge_keys[i]
        
        # Handle different merge key scenarios
        if isinstance(left_key, list) and isinstance(right_key, list):
            if len(left_key) == 1 and len(right_key) == 1:
                # If we have single keys, use them directly
                left_on = left_key[0]
                right_on = right_key[0]
                logger.info(f"Using left_on={left_on}, right_on={right_on} for merge")
                
                # Perform the merge
                result = result.merge(
                    df,
                    left_on=left_on,
                    right_on=right_on,
                    how=how,
                    suffixes=(suffixes[i-1], suffixes[i])
                )
                
                # Handle key column renaming
                if original_keys[i-1] == original_keys[i]:
                    # If both keys are the same, rename suffixed columns back to original
                    if f"{original_keys[i-1]}{suffixes[i-1]}" in result.columns:
                        result = result.rename(columns={f"{original_keys[i-1]}{suffixes[i-1]}": original_keys[i-1]})
                    if f"{original_keys[i]}{suffixes[i]}" in result.columns:
                        result = result.rename(columns={f"{original_keys[i]}{suffixes[i]}": original_keys[i]})
                else:
                    # If keys are different, keep both original names
                    if f"{left_on}{suffixes[i-1]}" in result.columns:
                        result = result.rename(columns={f"{left_on}{suffixes[i-1]}": left_on})
                    if f"{right_on}{suffixes[i]}" in result.columns:
                        result = result.rename(columns={f"{right_on}{suffixes[i]}": right_on})
            else:
                # If we have multiple keys, use them as lists
                logger.info(f"Using left_on={left_key}, right_on={right_key} for merge")
                result = result.merge(
                    df,
                    left_on=left_key,
                    right_on=right_key,
                    how=how,
                    suffixes=(suffixes[i-1], suffixes[i])
                )
        else:
            # If keys are strings or None, use them directly
            logger.info(f"Using on={left_key} for merge")
            result = result.merge(
                df,
                on=left_key,
                how=how,
                suffixes=(suffixes[i-1], suffixes[i])
            )
        logger.info(f"Shape after merge: {result.shape}")
    
    return result

def generate_metrics(df: pd.DataFrame, input_files: List[str], params: dict) -> dict:
    """Generate metrics about the merge operation.
    
    Args:
        df: Merged DataFrame
        input_files: List of input file paths
        params: Dictionary of merge parameters
        
    Returns:
        Dictionary of metrics
    """
    metrics = {
        'input_files': input_files,
        'parameters': params,
        'output_shape': {
            'rows': df.shape[0],
            'columns': df.shape[1]
        },
        'column_names': list(df.columns),
        'memory_usage_bytes': df.memory_usage(deep=True).sum()
    }
    
    return metrics

def main():
    # Get input files and output file
    input_files = snakemake.input.csv_files
    output_file = snakemake.output[0]
    metrics_file = snakemake.output.get('metrics', None)  # Get optional metrics file
    
    # Get parameters with defaults
    on = snakemake.params.get("on", None)
    how = snakemake.params.get("how", "inner")
    suffixes = snakemake.params.get("suffixes", None)
    
    params = {
        'on': on,
        'how': how,
        'suffixes': suffixes
    }
    
    try:
        # Perform the merge
        merged_df = merge_csv_files(
            input_files,
            on=on,
            how=how,
            suffixes=suffixes
        )
        
        # Save merged result
        merged_df.to_csv(output_file, index=False)
        logger.info(f"Saved merged result to {output_file}")
        
        # Generate and save metrics if requested
        if metrics_file:
            metrics = generate_metrics(merged_df, input_files, params)
            
            with open(metrics_file, 'w') as f:
                f.write("# CSV Merge Metrics\n")
                yaml.dump(metrics, f, default_flow_style=False)
            logger.info(f"Saved metrics to {metrics_file}")
            
    except Exception as e:
        logger.error(f"Error merging files: {e}")
        raise

if __name__ == "__main__":
    main() 