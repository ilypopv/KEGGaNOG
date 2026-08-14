"""Multi-sample three-panel heatmap generator module for KEGGaNOG profiles.

This module processes comprehensive matrix datasets across multiple sample columns,
dynamically computes proportional layout limits, and draws parallel subplots.
"""

import logging
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

_log = logging.getLogger(__name__)


def generate_heatmap_multi(
    kegg_decoder_file: pd.DataFrame,
    output_folder: str,
    dpi: int,
    color: str,
    figsize: tuple[float, float] | None = None,
) -> tuple[plt.Figure, Sequence[plt.Axes]]:
    """Generate a publication-grade dynamic wide three-panel heatmap for multiple samples.

    Processes pre-loaded multi-column functional matrices, computes scaled width ratios,
    splits pathway collections evenly, and renders integrated side-by-side matrices.

    Args:
        kegg_decoder_file: Input DataFrame containing multi-sample functional profiles.
        output_folder: Target location identifying active processing directory loops.
        dpi: Target resolution scale bounding the output drawing canvas.
        color: Target string colormap descriptor passed to downstream styling engines.
        figsize: Geometric allocation limits (width, height) defining canvas borders.

    Returns:
        tuple: Active Matplotlib Figure instance and the coordinate subplots Axes sequence.
    """
    _log.info("Generating heatmap...")

    # Validate sample elements and map numeric matrices under tracked progress
    with tqdm(total=3, desc="Preparing heatmap data") as pbar:
        df = pd.DataFrame(kegg_decoder_file)
        pbar.update(1)

        # Reshape layout spreadsheet structure via split operations
        df1, df2, df3 = split_dataframe_into_three_row_segments(df)
        pbar.update(2)

    # Establish dynamic proportional structural canvas limits if missing
    fig_w = 0.5 * (df.shape[1] - 2) + 19.5
    target_figsize = figsize if figsize is not None else (fig_w, 20.0)

    # Initialize structural subplots container canvas layers
    fig, axes, cbar_ax = create_three_panel_heatmap_figure(target_figsize)

    # Map underlying statistical matrices chunk-by-chunk onto partitioned views
    with tqdm(total=3, desc="Creating heatmap parts") as pbar:
        sns.heatmap(
            df1.set_index("Function"),
            cmap=color,
            annot=False,
            linewidths=0.5,
            ax=axes[0],
            cbar=False,
        )
        axes[0].set_title("Part 1")
        axes[0].tick_params(axis="x", rotation=45)
        axes[0].set_xticklabels(axes[0].get_xticklabels(), ha="right")
        pbar.update(1)

        sns.heatmap(
            df2.set_index("Function"),
            cmap=color,
            annot=False,
            linewidths=0.5,
            ax=axes[1],
            cbar=False,
        )
        axes[1].set_title("Part 2")
        axes[1].tick_params(axis="x", rotation=45)
        axes[1].set_xticklabels(axes[1].get_xticklabels(), ha="right")
        pbar.update(1)

        # Append shared structural vertical colorbars on the final pane
        sns.heatmap(
            df3.set_index("Function"),
            cmap=color,
            annot=False,
            linewidths=0.5,
            ax=axes[2],
            cbar_ax=cbar_ax,
            cbar_kws={"label": "Pathway completeness"},
        )
        axes[2].set_title("Part 3")
        axes[2].tick_params(axis="x", rotation=45)
        axes[2].set_xticklabels(axes[2].get_xticklabels(), ha="right")
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
