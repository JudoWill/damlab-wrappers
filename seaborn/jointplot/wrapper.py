"""Wrapper for seaborn jointplot functionality"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2025"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import yaml
import logging
from pathlib import Path
from typing import Optional, Dict, Any

# Configure logging to use Snakemake's log file
if "snakemake" in locals():
    log_file = snakemake.log[0] if hasattr(snakemake, 'log') and snakemake.log else None
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file
    )
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

logger = logging.getLogger('seaborn-jointplot-wrapper')

if "snakemake" not in locals():
    import snakemake  # type: ignore

def create_jointplot(df: pd.DataFrame, 
                    x: str,
                    y: str,
                    plot_kwargs: Optional[Dict[str, Any]] = None,
                    joint_kwargs: Optional[Dict[str, Any]] = None,
                    marginal_kwargs: Optional[Dict[str, Any]] = None,
                    fig_kwargs: Optional[Dict[str, Any]] = None) -> None:
    """Create a seaborn jointplot from a DataFrame.
    
    Args:
        df: Input DataFrame
        x: Column name for x-axis
        y: Column name for y-axis
        joint_kwargs: Arguments to pass to the joint plot
        marginal_kwargs: Arguments to pass to the marginal plots
        fig_kwargs: Arguments to pass to the figure
    """
    logger.info(f"Creating jointplot for {x} vs {y}")
    
    # Set up default parameters
    plot_kwargs = plot_kwargs or {}
    joint_kwargs = joint_kwargs or {}
    marginal_kwargs = marginal_kwargs or {}
    fig_kwargs = fig_kwargs or {}
    
    # Create the jointplot
    g = sns.jointplot(
        data=df,
        x=x,
        y=y,
        **plot_kwargs
    )
    
    # Apply marginal plot customizations if provided
    if marginal_kwargs:
        func_name = marginal_kwargs.pop('func', 'kdeplot')
        func = getattr(sns, func_name)
        g.plot_marginals(func, **marginal_kwargs)
    
    if joint_kwargs:
        func_name = joint_kwargs.pop('func', 'kdeplot')
        func = getattr(sns, func_name)
        g.plot_joint(func, **joint_kwargs)

    # Apply figure customizations if provided
    if fig_kwargs:
        g.figure.set(**fig_kwargs)
    
    # Save the figure
    g.savefig(snakemake.output[0])
    logger.info(f"Saved jointplot to {snakemake.output[0]}")
    
    # Close the figure to free memory
    plt.close()

def generate_metrics(df: pd.DataFrame, params: dict) -> dict:
    """Generate metrics about the jointplot operation.
    
    Args:
        df: Input DataFrame
        params: Dictionary of plot parameters
        
    Returns:
        Dictionary of metrics
    """
    metrics = {
        'input_shape': {
            'rows': df.shape[0],
            'columns': df.shape[1]
        },
        'column_names': list(df.columns),
        'parameters': params,
        'memory_usage_bytes': df.memory_usage(deep=True).sum()
    }
    
    return metrics

def main():
    # Get input file and output file
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    metrics_file = snakemake.output.get('metrics', None)  # Get optional metrics file
    
    # Get parameters
    x = snakemake.params['x']
    y = snakemake.params['y']
    plot_kwargs = snakemake.params.get('plot_kwargs', {})
    joint_kwargs = snakemake.params.get('joint_kwargs', {})
    marginal_kwargs = snakemake.params.get('marginal_kwargs', {})
    fig_kwargs = snakemake.params.get('fig_kwargs', {})
    fillna_kwargs = snakemake.params.get('fillna', None)  # Get optional fillna parameters
    
    params = {
        'x': x,
        'y': y,
        'joint_kwargs': joint_kwargs,
        'marginal_kwargs': marginal_kwargs,
        'fig_kwargs': fig_kwargs,
        'fillna': fillna_kwargs
    }
    
    try:
        # Read the input CSV file   
        logger.info(f"Reading data from {input_file}")
        
        df = pd.read_csv(input_file)
        logger.info(f"Read data with shape {df.shape}")

        # Apply fillna if specified
        if fillna_kwargs is not None:
            logger.info(f"Filling NA values with parameters: {fillna_kwargs}")
            df.fillna(inplace=True, **fillna_kwargs)
            logger.info("NA values filled")

        # Create the jointplot
        create_jointplot(
            df,
            x=x,
            y=y,
            plot_kwargs=plot_kwargs,
            joint_kwargs=joint_kwargs,
            marginal_kwargs=marginal_kwargs,
            fig_kwargs=fig_kwargs
        )
        
        # Generate and save metrics if requested
        if metrics_file:
            metrics = generate_metrics(df, params)
            
            with open(metrics_file, 'w') as f:
                f.write("# Seaborn Jointplot Metrics\n")
                yaml.dump(metrics, f, default_flow_style=False)
            logger.info(f"Saved metrics to {metrics_file}")
            
    except Exception as e:
        logger.error(f"Error creating jointplot: {e}")
        raise

if __name__ == "__main__":
    main() 