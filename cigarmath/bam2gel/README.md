# BAM2GEL

A wrapper for creating gel-like visualizations from BAM files using cigarmath. This tool processes one or more BAM files and creates a stacked histogram visualization that resembles an agarose gel, with each BAM file represented as a different lane.

## Features

- Process multiple BAM files in parallel
- Three different analysis modes:
  - `read_length`: Visualize read length distribution
  - `reference_block_size`: Visualize reference block size distribution
  - `query_block_size`: Visualize query block size distribution
- Memory-efficient streaming operations for large BAM files
- Customizable lane names and visualization parameters
- Outputs both visualization and detailed metrics

## Parameters

- `mode`: Analysis mode (default: "read_length")
  - Options: "read_length", "reference_block_size", "query_block_size"
- `names`: List of names for each lane (default: Sample_1, Sample_2, etc.)
- `bin_width`: Width of histogram bins in base pairs (default: 10)
- `max_size`: Maximum size to include in visualization (default: 1000)
- `min_size`: Minimum size to include in visualization (default: 1)
- `colorbar_min`: Minimum value for the color bar scale (default: 1)
- `colorbar_max`: Maximum value for the color bar scale (default: None, will use maximum observed value)
- `log_scale`: Whether to use logarithmic scale for the color bar (default: True)
- `fig_kwargs`: Dictionary of parameters passed to `plt.figure()` (default: {"figsize": (10, 6)})
- `plot_kwargs`: Dictionary of parameters passed to `sns.histplot()` (default: see below)
- `axis_kwargs`: Dictionary of parameters for axis customization (default: see below)

Default plot parameters:
```python
plot_kwargs = {
    "x": "lane",
    "y": "size",
    "bins": (len(names), np.arange(min_size, max_size + bin_width, bin_width)),
    "cbar": True,
    "cmap": "rocket_r",
    "vmin": colorbar_min,
    "vmax": colorbar_max
}
```

Default axis parameters:
```python
axis_kwargs = {
    "xlabel": "Lane",
    "ylabel": "Size (bp)",
    "title": f"Gel Visualization - {mode.replace('_', ' ').title()}"
}
```

## Outputs

1. PNG file containing the gel visualization
2. YAML file containing detailed metrics for each lane:
   - Total reads
   - Mean size
   - Median size
   - Minimum size
   - Maximum size

## Example Usage

Basic usage:
```python
rule bam2gel:
    input:
        bams=["sample1.bam", "sample2.bam"]
    output:
        png="gel_visualization.png",
        yaml="gel_metrics.yaml"
    params:
        mode="read_length",
        names=["Control", "Treatment"],
        bin_width=20,
        max_size=2000,
        min_size=1,
        colorbar_min=1,
        colorbar_max=1000,
        log_scale=True
    wrapper:
        "damlab-wrappers/cigarmath/bam2gel"
```

Customized visualization:
```python
rule bam2gel_custom:
    input:
        bams=["sample1.bam", "sample2.bam"]
    output:
        png="gel_visualization_custom.png",
        yaml="gel_metrics_custom.yaml"
    params:
        mode="read_length",
        names=["Control", "Treatment"],
        bin_width=20,
        max_size=2000,
        min_size=1,
        colorbar_min=1,
        colorbar_max=1000,
        log_scale=True,
        # Customize figure size and DPI
        fig_kwargs={
            "figsize": (12, 8),
            "dpi": 300
        },
        # Customize plot appearance
        plot_kwargs={
            "cmap": "viridis",  # Use a different colormap
            "cbar_kws": {"label": "Read Count"}  # Customize colorbar label
        },
        # Customize axis labels and title
        axis_kwargs={
            "xlabel": "Sample",
            "ylabel": "Fragment Length (bp)",
            "title": "Custom Gel Visualization"
        }
    wrapper:
        "damlab-wrappers/cigarmath/bam2gel"
```

## Dependencies

- cigarmath
- pandas
- seaborn
- matplotlib
- numpy
- pyyaml 