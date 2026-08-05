"""Single-sample three-panel heatmap generator module for KEGGaNOG profiles.

This module processes raw KEGG-Decoder output strings, normalizes single-vector
functional columns, and projects values onto a partitioned matrix grid layout.
"""

import warnings
from collections.abc import Sequence

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from .heatmaps_common import (
    create_three_panel_heatmap_figure,
    save_heatmap_png,
    split_dataframe_into_three_row_segments,
)


def generate_heatmap(
    kegg_decoder_file: str,
    output_folder: str,
    dpi: int,
    color: str,
    sample_name: str | None,
    figsize: tuple[float, float] | None = None,
    annot: bool = True,
) -> tuple[plt.Figure, Sequence[plt.Axes]]:
    """Generate a publication-grade partitioned three-panel heatmap for a single sample vector.

    Parses raw tabular streams, extracts targeted function column headers, maps
    completeness ranges, splits rows evenly, and paints unified side-by-side matrices.

    Args:
        kegg_decoder_file: System disk location pointing to raw text annotation matrices.
        output_folder: Target location identifying active processing directory loops.
        dpi: Target resolution scale bounding the output drawing canvas.
        color: Target string colormap descriptor passed to downstream styling engines.
        sample_name: Target sample identifier string mapping columns within data arrays.
        figsize: Geometric allocation limits (width, height) defining canvas borders.
        annot: Flag indicating whether cellular scalar elements are printed as text labels.

    Returns:
        tuple: Active Matplotlib Figure instance and the coordinate subplots Axes sequence.
    """
    # Parse and serialize disk spreadsheet rows into continuous strings
    with open(kegg_decoder_file, "r") as file:
        lines = file.readlines()

    # Validate sample elements and map numeric matrices under tracked progress
    with tqdm(total=3, desc="Preparing heatmap data") as pbar:
        header = lines[0].strip().split("\t")
        values = lines[1].strip().split("\t")

        # Fallback to general moniker if sample name context is absent
        target_column = sample_name if sample_name is not None else "Sample"

        data = {"Function": header[1:], target_column: [float(v) for v in values[1:]]}
        df = pd.DataFrame(data)
        pbar.update(1)

        # Reshape layout spreadsheet structure via split operations
        df1, df2, df3 = split_dataframe_into_three_row_segments(df)
        pbar.update(1)

        pbar.update(1)

    # Establish default structural canvas limits if missing
    target_figsize = figsize if figsize is not None else (20.0, 20.0)

    # Initialize structural subplots container canvas layers
    fig, axes, cbar_ax = create_three_panel_heatmap_figure(target_figsize)

    # Map underlying statistical matrices chunk-by-chunk onto partitioned views
    with tqdm(total=3, desc="Creating heatmap parts") as pbar:
        sns.heatmap(
            df1.pivot_table(values=target_column, index="Function", fill_value=0),
            cmap=color,
            annot=annot,
            linewidths=0.5,
            ax=axes[0],
            cbar=False,
        )
        axes[0].set_title("Part 1")
        pbar.update(1)

        sns.heatmap(
            df2.pivot_table(values=target_column, index="Function", fill_value=0),
            cmap=color,
            annot=annot,
            linewidths=0.5,
            ax=axes[1],
            cbar=False,
        )
        axes[1].set_title("Part 2")
        pbar.update(1)

        # Append shared structural vertical colorbars on the final pane
        sns.heatmap(
            df3.pivot_table(values=target_column, index="Function", fill_value=0),
            cmap=color,
            annot=annot,
            linewidths=0.5,
            ax=axes[2],
            cbar_ax=cbar_ax,
            cbar_kws={"label": "Pathway completeness"},
        )
        axes[2].set_title("Part 3")
        pbar.update(1)

        # Apply customized typography parameters to coordinate boundaries
        axes[1].set_ylabel("")
        axes[2].set_ylabel("")

    # Capture transient layout warnings and write final PNG assets to system disk
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=UserWarning, message=".*tight_layout.*"
        )
        plt.tight_layout(rect=(0.0, 0.0, 0.9, 1.0))

    save_heatmap_png(output_folder, dpi)

    # Return active layout context wrappers for external saving pipelines
    return fig, axes
