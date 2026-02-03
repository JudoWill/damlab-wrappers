"""Wrapper for seaborn plot functionality"""

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
from typing import Optional, Tuple, List, Union, Dict, Any

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

logger = logging.getLogger('seaborn-plot-wrapper')

if "snakemake" not in locals():
    import snakemake  # type: ignore

def create_seaborn_plot(df: pd.DataFrame, 
                       plot_type: Union[str, List[str]], 
                       plot_kwargs: Union[Dict[str, Any], List[Dict[str, Any]]],
                       fig_kwargs: Optional[Dict[str, Any]] = None,
                       axis_kwargs: Optional[Dict[str, Any]] = None) -> None:
    """Create one or more seaborn plots from a CSV file.
    
    Args:
        df: Input DataFrame
        plot_type: Type of seaborn plot to create (e.g., 'barplot', 'histplot', etc.)
                  or list of plot types for multiple plots
        plot_kwargs: Arguments to pass to the seaborn plot function(s)
        fig_kwargs: Arguments to pass to plt.figure()
        axis_kwargs: Arguments to pass to axis customization
    """
    logger.info(f"Creating plot(s) from {df.shape}")
    
    # Set up figure with custom parameters if provided
    fig_kwargs = fig_kwargs or {}
    fig, ax = plt.subplots(1,1, **fig_kwargs)
    
    # Convert single plot to list for consistent handling
    if isinstance(plot_type, str):
        plot_type = [plot_type]
    if isinstance(plot_kwargs, dict):
        plot_kwargs = [plot_kwargs]
    
    # Ensure plot_type and plot_kwargs lists have the same length
    if len(plot_type) != len(plot_kwargs):
        raise ValueError(f"Number of plot types ({len(plot_type)}) must match number of plot_kwargs ({len(plot_kwargs)})")
    
    # Create each plot
    for i, (ptype, pkwargs) in enumerate(zip(plot_type, plot_kwargs)):
        logger.info(f"Creating {ptype} plot {i+1}/{len(plot_type)}")
        
        # Get the plot function from seaborn
        plot_func = getattr(sns, ptype)
        
        # Create the plot
        plot_func(data=df, ax=ax, **pkwargs)
    
    # Apply axis customizations if provided
    axis_kwargs = axis_kwargs or {}
    
    # Handle xlabel
    if 'xlabel' in axis_kwargs:
        ax.set_xlabel(axis_kwargs.pop('xlabel'))
    
    # Handle ylabel
    if 'ylabel' in axis_kwargs:
        ax.set_ylabel(axis_kwargs.pop('ylabel'))
    
    # Handle xlim
    if 'xlim' in axis_kwargs:
        ax.set_xlim(axis_kwargs.pop('xlim'))
    
    # Handle ylim
    if 'ylim' in axis_kwargs:
        ax.set_ylim(axis_kwargs.pop('ylim'))
    
    # Handle despine
    if 'despine' in axis_kwargs:
        despine_params = axis_kwargs.pop('despine')
        if isinstance(despine_params, bool):
            sns.despine() if despine_params else None
        else:
            sns.despine(**despine_params)
    
    # Handle legend
    if 'legend' in axis_kwargs:
        legend_param = axis_kwargs.pop('legend')
        if isinstance(legend_param, bool):
            if legend_param:
                ax.legend()
            else:
                ax.legend().remove()
        else:
            ax.legend(**legend_param)
    
    # Apply any remaining axis customizations
    for key, value in axis_kwargs.items():
        if hasattr(ax, f'set_{key}'):
            getattr(ax, f'set_{key}')(value)
    
    # Save the figure
    fig.savefig(snakemake.output[0])
    logger.info(f"Saved plot to {snakemake.output[0]}")
    
    # Close the figure to free memory
    plt.close()

def generate_metrics(df: pd.DataFrame, params: dict) -> dict:
    """Generate metrics about the plot operation.
    
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
    plot_type = snakemake.params['plot']
    plot_kwargs = snakemake.params.get('plot_kwargs', {})
    fig_kwargs = snakemake.params.get('fig_kwargs', {})
    axis_kwargs = snakemake.params.get('axis_kwargs', {})
    fillna_kwargs = snakemake.params.get('fillna', None)  # Get optional fillna parameters
    
    # Convert plot_type and plot_kwargs to lists if they aren't already
    if not isinstance(plot_type, list):
        plot_type = [plot_type]
    if not isinstance(plot_kwargs, list):
        plot_kwargs = [plot_kwargs]
    
    # Ensure plot_type and plot_kwargs lists have the same length
    if len(plot_type) != len(plot_kwargs):
        raise ValueError(f"Number of plot types ({len(plot_type)}) must match number of plot_kwargs ({len(plot_kwargs)})")
    
    params = {
        'plot_type': plot_type,
        'plot_kwargs': plot_kwargs,
        'fig_kwargs': fig_kwargs,
        'axis_kwargs': axis_kwargs,
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

        # Create the plot
        create_seaborn_plot(
            df,
            plot_type=plot_type,
            plot_kwargs=plot_kwargs,
            fig_kwargs=fig_kwargs,
            axis_kwargs=axis_kwargs
        )
        
        # Generate and save metrics if requested
        if metrics_file:
            metrics = generate_metrics(df, params)
            
            with open(metrics_file, 'w') as f:
                f.write("# Seaborn Plot Metrics\n")
                yaml.dump(metrics, f, default_flow_style=False)
            logger.info(f"Saved metrics to {metrics_file}")
            
    except Exception as e:
        logger.error(f"Error creating plot: {e}")
        raise

if __name__ == "__main__":
    main()
