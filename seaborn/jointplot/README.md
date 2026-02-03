# Seaborn Jointplot Wrapper

A Snakemake wrapper for creating seaborn jointplots from CSV data.

## Overview

This wrapper provides a simple interface to create seaborn jointplots from CSV data. Jointplots combine scatter plots with histograms or kernel density estimates on the margins, providing a comprehensive view of the relationship between two variables. The wrapper supports all jointplot parameters and allows customization of both the main plot and marginal plots. It also supports filling NA values in the data before plotting.

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
rule create_jointplot:
    input:
        "data.csv"
    output:
        "jointplot.png",
        metrics="jointplot_metrics.yaml"
    params:
        x="x_column",  # Column name for x-axis
        y="y_column",  # Column name for y-axis
        joint_kwargs={  # Arguments to pass to the joint plot
            "kind": "scatter",
            "hue": "group",
            "alpha": 0.5
        },
        marginal_kwargs={  # Arguments to pass to the marginal plots
            "kde": True,
            "bins": 30
        },
        fig_kwargs={  # Arguments to pass to the figure
            "figsize": (10, 10),
            "dpi": 300
        },
        fillna={  # Optional arguments for pandas fillna
            "value": 0,  # Fill NA with 0
            "method": "ffill"  # Forward fill
        }
    wrapper:
        "file:damlab-wrappers/seaborn/jointplot"
```

### Parameters

- `x`: Column name for x-axis
- `y`: Column name for y-axis
- `joint_kwargs`: Dictionary of arguments to pass to the joint plot:
  - `kind`: Type of plot to draw ('scatter', 'kde', 'hist', 'hex', 'reg', 'resid')
  - `hue`: Variable that determines the color of plot elements
  - `alpha`: Transparency level for points
  - Any other valid seaborn.jointplot() parameters
- `marginal_kwargs`: Dictionary of arguments to pass to the marginal plots:
  - `kde`: Whether to plot a kernel density estimate
  - `hist`: Whether to plot a histogram
  - `bins`: Number of bins for histogram
  - Any other valid seaborn.kdeplot() or seaborn.histplot() parameters
- `fig_kwargs`: Dictionary of arguments to pass to the figure:
  - `figsize`: Tuple of (width, height) in inches
  - `dpi`: Dots per inch for the figure
  - Any other valid plt.figure() parameters
- `fillna`: Optional dictionary of arguments to pass to pandas.DataFrame.fillna():
  - `value`: Scalar value to use to fill NA values
  - `method`: Method to use for filling holes in reindexed Series ('ffill', 'bfill', etc.)
  - `axis`: Axis along which to fill missing values
  - `inplace`: Whether to modify the DataFrame in place (always True in this wrapper)
  - Any other valid pandas.DataFrame.fillna() parameters

### Outputs

- The main output is the jointplot file (typically a PNG or PDF)
- An optional metrics file (YAML format) containing information about the input data and plot parameters

## Examples

### Basic Jointplot

```python
rule create_basic_jointplot:
    input:
        "data.csv"
    output:
        "jointplot.png"
    params:
        x="x_values",
        y="y_values",
        joint_kwargs={
            "kind": "scatter",
            "alpha": 0.5
        }
    wrapper:
        "file:damlab-wrappers/seaborn/jointplot"
```

### Jointplot with KDE Marginals

```python
rule create_kde_jointplot:
    input:
        "data.csv"
    output:
        "kde_jointplot.png"
    params:
        x="x_values",
        y="y_values",
        joint_kwargs={
            "kind": "kde",
            "hue": "group"
        },
        marginal_kwargs={
            "kde": True,
            "fill": True
        }
    wrapper:
        "file:damlab-wrappers/seaborn/jointplot"
```

### Jointplot with NA Handling

```python
rule create_jointplot_with_na:
    input:
        "data.csv"
    output:
        "jointplot.png"
    params:
        x="x_values",
        y="y_values",
        joint_kwargs={
            "kind": "scatter",
            "hue": "group"
        },
        fillna={
            "value": 0,  # Fill NA with 0
            "method": "ffill"  # Forward fill remaining NA
        }
    wrapper:
        "file:damlab-wrappers/seaborn/jointplot"
```

## Version History

- 1.0.0: Initial release 