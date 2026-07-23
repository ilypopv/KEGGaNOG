#!/usr/bin/env python3
"""Multi-sample grouped three-panel heatmap generator module for KEGGaNOG.

This module processes composite matrices tracking functional group buckets across multiple
samples, calculates custom geometric margin pads, and draws unified partitioned panels.
"""

import warnings
from typing import Dict, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from .grouped_heatmap import function_groups
from .heatmaps_common import (
    GROUPED_PART1_GROUPS,
    GROUPED_PART2_GROUPS,
    GROUPED_PART3_GROUPS,
    insert_split_rows_between_groups,
    save_heatmap_png,
)


def generate_grouped_heatmap_multi(
    kegg_decoder_file: pd.DataFrame,
    output_folder: str,
    dpi: int,
    color: str,
    figsize: Optional[Tuple[float, float]] = None,
) -> Tuple[plt.Figure, Sequence[plt.Axes]]:
    """Generate a functional grouped dynamic wide three-panel heatmap for multiple samples.

    Processes multidimensional pathway matrices, dynamically evaluates category margins,
    injects row spacers, draws parallel subplots, and formats metadata bounding boxes.

    Args:
        kegg_decoder_file: Input DataFrame containing multi-sample functional profiles.
        output_folder: Target location identifying active processing directory loops.
        dpi: Target resolution scale bounding the output drawing canvas.
        color: Target string colormap descriptor passed to downstream styling engines.
        figsize: Geometric allocation limits (width, height) defining canvas borders.

    Returns:
        tuple: Active Matplotlib Figure instance and the coordinate subplots Axes sequence.
    """
    # Validate sample elements and map numeric matrices under tracked progress
    with tqdm(total=6, desc="Preparing heatmap data") as pbar:
        function_groups_lower = {
            group: {func.lower() for func in funcs}
            for group, funcs in function_groups.items()
        }

        kegg_decoder_file["Group"] = kegg_decoder_file["Function"].apply(
            lambda x: next(
                (
                    group
                    for group, funcs in function_groups_lower.items()
                    if x.lower() in funcs
                ),
                "Miscellaneous",
            )
        )
        pbar.update(1)

        kegg_decoder_file = kegg_decoder_file.sort_values(
            by=["Group", "Function"]
        ).reset_index(drop=True)
        pbar.update(1)

        kegg_decoder_file["Function"] = pd.Categorical(
            kegg_decoder_file["Function"],
            categories=kegg_decoder_file["Function"],
            ordered=True,
        )
        pbar.update(1)

        part1_groups = GROUPED_PART1_GROUPS
        part2_groups = GROUPED_PART2_GROUPS
        part3_groups = GROUPED_PART3_GROUPS

        # Reshape layout spreadsheet structure via split operations
        part1 = kegg_decoder_file[
            kegg_decoder_file["Group"].isin(part1_groups)
        ].reset_index(drop=True)
        pbar.update(1)
        part2 = kegg_decoder_file[
            kegg_decoder_file["Group"].isin(part2_groups)
        ].reset_index(drop=True)
        pbar.update(1)
        part3 = kegg_decoder_file[
            kegg_decoder_file["Group"].isin(part3_groups)
        ].reset_index(drop=True)
        pbar.update(1)

    # Insert artificial NaN spacer rows between category groups
    with tqdm(total=6, desc="Adding split between groups") as pbar:
        part1 = insert_split_rows_between_groups(
            kegg_decoder_file[kegg_decoder_file["Group"].isin(part1_groups)],
            part1_groups,
        ).reset_index(drop=True)
        pbar.update(1)
        part2 = insert_split_rows_between_groups(
            kegg_decoder_file[kegg_decoder_file["Group"].isin(part2_groups)],
            part2_groups,
        ).reset_index(drop=True)
        pbar.update(1)
        part3 = insert_split_rows_between_groups(
            kegg_decoder_file[kegg_decoder_file["Group"].isin(part3_groups)],
            part3_groups,
        ).reset_index(drop=True)
        pbar.update(1)

        part1["Function"] = pd.Categorical(
            part1["Function"], categories=part1["Function"], ordered=True
        )
        pbar.update(1)
        part2["Function"] = pd.Categorical(
            part2["Function"], categories=part2["Function"], ordered=True
        )
        pbar.update(1)
        part3["Function"] = pd.Categorical(
            part3["Function"], categories=part3["Function"], ordered=True
        )
        pbar.update(1)

    # Establish dynamic proportional structural canvas limits if missing
    fig_w = 0.5 * (part1.shape[1] - 2) + 27.5
    target_figsize = figsize if figsize is not None else (fig_w, 20.0)

    # Mathematical polynomial margin calibration mapping input wide metrics safely
    left_pad = (
        0.20 - 0.0020833 * (part1.shape[1] - 2) + 0.0020833 * (part1.shape[1] - 2) ** 2
    )

    # Initialize structural subplots container canvas layers
    fig, axes_array = plt.subplots(1, 3, figsize=target_figsize)
    axes: Sequence[plt.Axes] = (
        axes_array.tolist() if hasattr(axes_array, "tolist") else axes_array
    )

    cbar_ax = fig.add_axes((0.92, 0.4, 0.02, 0.2))

    plt.subplots_adjust(left=left_pad, right=0.85, wspace=0.4)

    # Functional label generator context helpers
    def add_group_labels(
        current_ax: plt.Axes, part_df: pd.DataFrame, group_labels: Sequence[str]
    ) -> None:
        for group in group_labels:
            group_indices = np.where(part_df["Group"] == group)[0]
            if len(group_indices) > 0:
                y_position = float(np.mean(group_indices) + 0.5)
                x_position = -(left_pad / 2)
                current_ax.text(
                    x_position,
                    y_position,
                    group,
                    fontsize=12,
                    ha="right",
                    va="center",
                    weight="bold",
                    bbox=dict(
                        boxstyle="round,pad=0.3", edgecolor="none", facecolor="white"
                    ),
                )

    def plot_heatmap(
        part_df: pd.DataFrame,
        group_labels: Sequence[str],
        ax: plt.Axes,
        cbar: bool,
        cbar_axis: Optional[plt.Axes] = None,
        cbar_kws: Optional[Dict[str, str]] = None,
    ) -> None:
        value_columns = part_df.columns[1:-1]
        part_df[value_columns] = part_df[value_columns].fillna(0)

        pivot_table = part_df.set_index("Function")[value_columns]
        mask = pivot_table.index.str.startswith("split_")

        sns.heatmap(
            pivot_table,
            cmap=color,
            annot=False,
            linewidths=0.5,
            ax=ax,
            cbar=cbar,
            cbar_ax=cbar_axis,
            cbar_kws=cbar_kws,
            mask=np.tile(mask[:, None], (1, pivot_table.shape[1])),
        )
        # Derive tick positions directly from the index — get_yticklabels()
        # returns empty text before the figure is rendered, making the old
        # set_visible(False) loop a silent no-op.
        non_split_pos = [
            i + 0.5
            for i, name in enumerate(pivot_table.index)
            if not str(name).startswith("split_")
        ]
        non_split_labels = [
            name for name in pivot_table.index if not str(name).startswith("split_")
        ]
        ax.set_yticks(non_split_pos)
        ax.set_yticklabels(non_split_labels, rotation=0)
        ax.tick_params(axis="y", which="both", left=False)
        add_group_labels(ax, part_df, group_labels)

    # Map underlying statistical matrices chunk-by-chunk onto partitioned views
    with tqdm(total=3, desc="Creating heatmap parts") as pbar:
        plot_heatmap(part1, part1_groups, axes[0], cbar=False)
        axes[0].set_title("Part 1")
        axes[0].tick_params(axis="x", rotation=45)
        axes[0].set_xticklabels(axes[0].get_xticklabels(), ha="right")
        pbar.update(1)

        plot_heatmap(part2, part2_groups, axes[1], cbar=False)
        axes[1].set_title("Part 2")
        axes[1].tick_params(axis="x", rotation=45)
        axes[1].set_xticklabels(axes[1].get_xticklabels(), ha="right")
        pbar.update(1)

        plot_heatmap(
            part3,
            part3_groups,
            axes[2],
            cbar=True,
            cbar_axis=cbar_ax,
            cbar_kws={"label": "Pathway completeness"},
        )
        axes[2].set_title("Part 3")
        axes[2].tick_params(axis="x", rotation=45)
        axes[2].set_xticklabels(axes[2].get_xticklabels(), ha="right")
        pbar.update(1)

    # Apply customized typography parameters to coordinate boundaries
    axes[0].set_ylabel("")
    axes[1].set_ylabel("")
    axes[2].set_ylabel("")

    for ax in axes:
        ax.yaxis.tick_right()
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, va="center", ha="left")

    # Apply transient layout warnings and write final PNG assets to system disk
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=UserWarning, message=".*tight_layout.*"
        )
        plt.tight_layout(rect=(0.0, 0.0, 0.9, 1.0))

    save_heatmap_png(output_folder, dpi)

    # Return active layout context wrappers for external saving pipelines
    return fig, axes
