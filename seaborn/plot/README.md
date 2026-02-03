# Seaborn Plot Wrapper

A Snakemake wrapper for creating seaborn plots from CSV data.

## Overview

This wrapper provides a simple interface to create various types of seaborn plots from CSV data. It supports all seaborn plot types (barplot, histplot, scatterplot, etc.) and allows passing any additional parameters to the underlying seaborn functions. It also provides control over figure and axis properties. Multiple plots can be created on the same figure by providing lists of plot types and parameters. The wrapper also supports filling NA values in the data before plotting.

## Installation

The wrapper requires the following dependencies:
- pandas
- seaborn
- matplotlib
- pyyaml

These can be installed using the provided environment.yaml file:

```bash
conda env create -f environment.yaml
```

## Usage

### In Snakemake

```python
rule create_plot:
    input:
        "data.csv"
    output:
        "plot.png",
        metrics="plot_metrics.yaml"
    params:
        plot="barplot",  # Type of seaborn plot to create (or list of plot types)
        plot_kwargs={    # Arguments to pass to the seaborn plot function (or list of kwargs)
            "x": "category",
            "y": "value",
            "hue": "group"
        },
        fig_kwargs={    # Arguments to pass to plt.figure()
            "figsize": (10, 6),
            "dpi": 300
        },
        axis_kwargs={   # Arguments for axis customization
            "xlabel": "Category",
            "ylabel": "Value",
            "xlim": [0, 10],
            "ylim": [0, 100],
            "despine": True,
            "legend": True
        },
        fillna={        # Optional arguments for pandas fillna
            "value": 0,  # Fill NA with 0
            "method": "ffill"  # Forward fill
        }
    wrapper:
        "file:damlab-wrappers/seaborn/plot"
```

### Parameters

- `plot`: The type of seaborn plot to create (e.g., 'barplot', 'histplot', 'scatterplot', etc.) or a list of plot types for multiple plots
- `plot_kwargs`: Dictionary of arguments to pass to the seaborn plot function or a list of dictionaries for multiple plots
- `fig_kwargs`: Dictionary of arguments to pass to plt.figure():
  - `figsize`: Tuple of (width, height) in inches
  - `dpi`: Dots per inch for the figure
  - Any other valid plt.figure() parameters
- `axis_kwargs`: Dictionary of arguments for axis customization:
  - `xlabel`: X-axis label
  - `ylabel`: Y-axis label
  - `xlim`: X-axis limits as [min, max]
  - `ylim`: Y-axis limits as [min, max]
  - `despine`: Boolean or dictionary of arguments for sns.despine()
  - `legend`: Boolean or dictionary of arguments for ax.legend()
  - Any other valid axis setter methods (e.g., 'title', 'grid', etc.)
- `fillna`: Optional dictionary of arguments to pass to pandas.DataFrame.fillna():
  - `value`: Scalar value to use to fill NA values
  - `method`: Method to use for filling holes in reindexed Series ('ffill', 'bfill', etc.)
  - `axis`: Axis along which to fill missing values
  - `inplace`: Whether to modify the DataFrame in place (always True in this wrapper)
  - Any other valid pandas.DataFrame.fillna() parameters

### Outputs

- The main output is the plot file (typically a PNG or PDF)
- An optional metrics file (YAML format) containing information about the input data and plot parameters

## Examples

### Basic Barplot

```python
rule create_barplot:
    input:
        "data.csv"
    output:
        "barplot.png"
    params:
        plot="barplot",
        plot_kwargs={
            "x": "category",
            "y": "value",
            "hue": "group",
            "ci": "sd"
        }
    wrapper:
        "file:damlab-wrappers/seaborn/plot"
```

### Plot with NA Handling

```python
rule create_plot_with_na:
    input:
        "data.csv"
    output:
        "plot.png"
    params:
        plot="barplot",
        plot_kwargs={
            "x": "category",
            "y": "value",
            "hue": "group"
        },
        fillna={
            "value": 0,  # Fill NA with 0
            "method": "ffill"  # Forward fill remaining NA
        }
    wrapper:
        "file:damlab-wrappers/seaborn/plot"
```

### Multiple Plots on Same Figure

```python
rule create_complex_plot:
    input:
        "data.csv"
    output:
        "complex_plot.png"
    params:
        plot=["scatterplot", "lineplot"],  # Create both a scatter and line plot
        plot_kwargs=[
            {  # Arguments for scatter plot
                "x": "x_values",
                "y": "y_values",
                "hue": "group",
                "alpha": 0.5
            },
            {  # Arguments for line plot
                "x": "x_values",
                "y": "y_values",
                "hue": "group",
                "style": "group",
                "markers": True
            }
        ],
        fig_kwargs={
            "figsize": (12, 8),
            "dpi": 300
        },
        axis_kwargs={
            "xlabel": "X Values",
            "ylabel": "Y Values",
            "title": "Combined Scatter and Line Plot",
            "despine": True,
            "legend": {"title": "Groups", "loc": "upper right"}
        },
        fillna={
            "value": 0
        }
    wrapper:
        "file:damlab-wrappers/seaborn/plot"
```

## Version History

- 1.0.0: Initial release 