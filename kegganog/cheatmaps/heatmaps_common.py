"""Shared helper utilities and structural constants for KEGGaNOG heatmap layouts.

This module houses categorical functional group buckets and geometric partitioning
algorithms utilizing multi-panel subplots to visualize KEGG module profiles.
"""

from __future__ import annotations

import os
from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

# Group-category buckets for three-panel grouped heatmaps (must stay in lockstep).
GROUPED_PART1_GROUPS: list[str] = [
    "Amino acid metabolism",
    "Arsenic reduction",
    "Bacterial secretion systems",
    "Biofilm formation",
    "Carbohydrate metabolism",
    "Photosynthesis",
]
GROUPED_PART2_GROUPS: list[str] = [
    "Carbon degradation",
    "Carbon fixation",
    "Cell mobility",
    "Genetic competence",
    "Hydrogen redox",
    "Metal transporters",
    "Methanogenesis",
    "Miscellaneous",
]
GROUPED_PART3_GROUPS: list[str] = [
    "Nitrogen metabolism",
    "Oxidative phosphorylation",
    "Sulfur metabolism",
    "Transporters",
    "Vitamin biosynthesis",
]


def split_dataframe_into_three_row_segments(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Split pathway records into three contiguous vertical blocks using floor division.

    Args:
        df: Input DataFrame containing formatted completeness vectors.

    Returns:
        tuple: Matrix chunks representing partitioned rows for multi-panel rendering.
    """
    # Compute matrix partition bounds based on absolute row volume
    num_rows = len(df)
    split_size = num_rows // 3

    # Slice target positions isolating contiguous horizontal slices
    return (
        df.iloc[:split_size],
        df.iloc[split_size : 2 * split_size],
        df.iloc[2 * split_size :],
    )


def insert_split_rows_between_groups(
    df: pd.DataFrame, groups: Sequence[str]
) -> pd.DataFrame:
    """Insert artificial NaN spacer rows between discrete functional categories.

    Args:
        df: Input DataFrame containing functional group reference indices.
        groups: Categorical collection sequence dictating the validation block order.

    Returns:
        pd.DataFrame: Augmented matrix layout containing explicit group boundaries.
    """
    # Iterate through categories and group target row data blocks
    new_rows: list[pd.DataFrame] = []
    for group in groups:
        group_rows = df[df["Group"] == group]
        new_rows.append(group_rows)

        # Conditionally append dummy empty rows to act as visual grid splitters
        if group != groups[-1]:
            empty_row = pd.DataFrame(
                [["split_" + f"{group}"] + [np.nan] * (df.shape[1] - 1)],
                columns=df.columns,
            )
            new_rows.append(empty_row)

    # Reconstruct continuous matrix framing optimized categories
    return pd.concat(new_rows, ignore_index=True)


def create_three_panel_heatmap_figure(
    figsize: tuple[float, float],
) -> tuple[plt.Figure, Sequence[plt.Axes], plt.Axes]:
    """Initialize a multi-panel grid layout consisting of three aligned subplots and a colorbar.

    Args:
        figsize: Geometric allocation limits (width, height) defining canvas borders.

    Returns:
        tuple: Matplotlib figure container, structured subplots axes collection, and colorbar axes.
    """
    import matplotlib.pyplot as plt

    # Initialize structural subplots container canvas layers
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # Manually position a fixed colorbar widgets coordinate layout frame
    cbar_ax = fig.add_axes((0.92, 0.4, 0.02, 0.2))

    return fig, axes, cbar_ax


def save_heatmap_png(output_folder: str, dpi: int) -> str:
    """Serialize the current drawing canvas into a static high-resolution PNG asset.

    Args:
        output_folder: Destination path identifying target workspace directories.
        dpi: Target pixel resolution scale bounding the raster layout grid.

    Returns:
        str: Absolute destination path mapping the exported PNG asset file.
    """
    import matplotlib.pyplot as plt

    # Build file output pathways and update transactional context tracking bars
    output_file = os.path.join(output_folder, "heatmap_figure.png")

    with tqdm(total=1, desc="Saving plot") as pbar:
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
        pbar.update(1)

    # Flush current drawing stream pipelines and return active path locations
    plt.show()
    return output_file
