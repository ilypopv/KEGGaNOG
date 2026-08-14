"""Streamgraph visualization module for KEGG module completion profiles.

This module builds continuous, multi-component streamgraph layouts tracking
aggregated functional pathway group completeness variations derived from
eggNOG-mapper orthology annotations across independent samples.
"""

from collections.abc import Sequence
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ..cheatmaps.grouped_heatmap import function_groups
from .base import KgnPlotBase


class KgnStreamgraph(KgnPlotBase):
    """Orchestration context wrapper encapsulating Matplotlib streamgraph layouts."""

    def __init__(self, fig: plt.Figure, ax: plt.Axes) -> None:
        """Initialize the streamgraph canvas with layout metrics.

        Args:
            fig: The Matplotlib Figure container hosting the drawing canvas.
            ax: The core underlying Axes coordinate grid mapper.
        """
        super().__init__(fig, ax)


def streamgraph(
    df: pd.DataFrame,
    figsize: tuple[float, float] = (14.0, 7.0),
    cmap: str | Sequence[str] = "tab20",
    bar_width: float = 0.6,
    fill_alpha: float = 1.0,
    edgecolor: str | None = None,
    edge_linewidth: float = 0.3,
    title: str | None = None,
    title_fontsize: float = 16.0,
    title_color: str = "black",
    title_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    title_style: Literal["normal", "italic", "oblique"] = "normal",
    xlabel: str = "Samples",
    xlabel_fontsize: float = 12.0,
    xlabel_color: str = "black",
    xlabel_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    xlabel_style: Literal["normal", "italic", "oblique"] = "normal",
    ylabel: str = "Total Completeness",
    ylabel_fontsize: float = 12.0,
    ylabel_color: str = "black",
    ylabel_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    ylabel_style: Literal["normal", "italic", "oblique"] = "normal",
    xticks_rotation: float = 0.0,
    xticks_ha: Literal["left", "right", "center"] = "center",
    xticks_fontsize: float = 12.0,
    xticks_color: str = "black",
    xticks_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    xticks_style: Literal["normal", "italic", "oblique"] = "normal",
    background_color: str | None = "white",
    grid: bool = True,
    grid_linestyle: str = "--",
    grid_alpha: float = 0.7,
    legend_fontsize: float = 9.0,
    legend_loc: str = "upper left",
    legend_bbox: tuple[float, float] = (1.05, 1.0),
    show_legend: bool = True,
) -> KgnStreamgraph:
    """Generate a publication-grade customizable streamgraph for KEGG pathway groups.

    Transforms eggNOG-inferred functional profiles into aligned cross-pivoted
    matrix structures, aggregates completeness scores based on predefined functional
    categories, and generates smooth stream/stackplot transitions across analyzed samples.

    Args:
        df: Input DataFrame containing mapped functional metrics ('Function' and sample columns).
        figsize: Geometric allocation limits (width, height) defining canvas borders.
        cmap: Target string name lookup or a sequential list of direct hexadecimal colors.
        bar_width: Horizontal width metric allocated to independent drawing column bars.
        fill_alpha: Transparency scale constraint applied to drawing layers.
        edgecolor: Border styling color outline separating adjacent stacked blocks.
        edge_linewidth: Thickness parameter of border outlines enclosing drawing cells.
        title: Global text message identifier rendering above the drawing matrix.
        title_fontsize, title_color, title_weight, title_style: Font properties for title.
        xlabel, ylabel: Text content values mapping coordinates descriptors.
        xlabel_fontsize, ylabel_fontsize: Typography scale indices.
        xlabel_color, ylabel_color: Text color variables mapping target labels.
        xlabel_weight, ylabel_weight: Structural typographic density metrics.
        xlabel_style, ylabel_style: Geometric font slope configurations.
        xticks_rotation, xticks_ha: Position variables mapping target X tick attributes.
        xticks_fontsize, xticks_color, xticks_weight, xticks_style: X-tick typography rules.
        background_color: Primary layout canvas backdrop color mapping.
        grid: Toggles background coordinate reference line structures.
        grid_linestyle: Grid line texture rendering parameter.
        grid_alpha: Opacity index managing visibility bounds of grid elements.
        legend_fontsize: Typography density constraints mapped into explicit sub-labels.
        legend_loc: Positional anchoring code identifier tracking layout widgets.
        legend_bbox: Coordinate anchor box offsets defining bounding regions for legends.
        show_legend: If False, completely suppresses widget layer execution.

    Returns:
        KgnStreamgraph: Container instance holding references to optimized figures.
    """
    function_to_group = {}
    for group, funcs in function_groups.items():
        for func in funcs:
            function_to_group[func.lower()] = group

    # Validate sample elements and apply strict ordered categorical indices
    working_df = df.copy()
    working_df["Group"] = working_df["Function"].str.lower().map(function_to_group)
    df_grouped = working_df.dropna(subset=["Group"])
    df_grouped_sum = df_grouped.groupby("Group").sum(numeric_only=True)

    # Reshape layout spreadsheet structure via pivot operations
    df_plot = df_grouped_sum.T

    # Establish palette map dictionaries compliant with static analysis
    if isinstance(cmap, str):
        colors = sns.color_palette(cmap, n_colors=len(df_plot.columns))
    elif isinstance(cmap, (list, tuple, np.ndarray, pd.Series)) or hasattr(
        cmap, "__iter__"
    ):
        colors = cmap
    else:
        colors = sns.color_palette("tab20", n_colors=len(df_plot.columns))

    # Initialize structural subplots container canvas layers
    fig, ax = plt.subplots(figsize=figsize, facecolor=background_color)

    centers = np.arange(len(df_plot.index))
    xs = np.column_stack((centers - bar_width / 2, centers + bar_width / 2)).flatten()
    ys = np.repeat(df_plot.values.T, 2, axis=1)

    layers = ax.stackplot(
        xs,
        ys,
        labels=df_plot.columns,
        colors=colors,
        alpha=fill_alpha,
        zorder=3,
    )

    if edgecolor is not None:
        for poly in layers:
            poly.set_edgecolor(edgecolor)
            poly.set_linewidth(edge_linewidth)

    # Apply customized typography parameters to coordinate boundaries
    if title:
        ax.set_title(
            title,
            fontsize=title_fontsize,
            color=title_color,
            weight=title_weight,
            style=title_style,
        )
    ax.set_xlabel(
        xlabel,
        fontsize=xlabel_fontsize,
        color=xlabel_color,
        weight=xlabel_weight,
        style=xlabel_style,
    )
    ax.set_ylabel(
        ylabel,
        fontsize=ylabel_fontsize,
        color=ylabel_color,
        weight=ylabel_weight,
        style=ylabel_style,
    )

    # Build legend overlays if requested by execution flags
    if show_legend:
        fig.subplots_adjust(right=0.72)
        ax.legend(
            title="Pathway Group",
            bbox_to_anchor=legend_bbox,
            loc=legend_loc,
            fontsize=legend_fontsize,
            borderaxespad=0.5,
        )

    if grid:
        ax.grid(axis="x", linestyle=grid_linestyle, alpha=grid_alpha, zorder=0)

    # Configure ticks geometric parameters and close active canvas stream descriptors
    plt.xticks(
        centers,
        df_plot.index.tolist(),
        rotation=xticks_rotation,
        ha=xticks_ha,
        fontsize=xticks_fontsize,
        color=xticks_color,
        weight=xticks_weight,
        style=xticks_style,
    )

    ax.set_xlim(-0.5, len(df_plot.index) - 0.5)

    plt.close(fig)

    return KgnStreamgraph(fig, ax)
