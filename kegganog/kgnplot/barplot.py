"""Horizontal barplot visualization module for KEGG module completion profiles.

This module builds standalone sorted horizontal bar charts tracking individual
pathway completeness scores extracted from eggNOG-mapper functional annotations.
"""

from collections.abc import Sequence
from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .base import KgnPlotBase


class KgnBarplot(KgnPlotBase):
    """Orchestration context wrapper encapsulating Matplotlib barplot layouts."""

    def __init__(self, fig: plt.Figure, ax: plt.Axes) -> None:
        """Initialize the barplot canvas with layout metrics.

        Args:
            fig: The Matplotlib Figure container hosting the drawing canvas.
            ax: The core underlying Axes coordinate grid mapper.
        """
        super().__init__(fig, ax)


def barplot(
    df: pd.DataFrame,
    figsize: tuple[float, float] = (8.0, 12.0),
    cmap: str | Sequence[str] = "Greens",
    cmap_range: tuple[int, int] = (8, 30),
    title: str | None = None,
    title_fontsize: float = 16.0,
    title_color: str = "black",
    title_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    title_style: Literal["normal", "italic", "oblique"] = "normal",
    xlabel: str = "Pathway completeness",
    xlabel_fontsize: float = 14.0,
    xlabel_color: str = "black",
    xlabel_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    xlabel_style: Literal["normal", "italic", "oblique"] = "normal",
    ylabel: str = "Pathway",
    ylabel_fontsize: float = 14.0,
    ylabel_color: str = "black",
    ylabel_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    ylabel_style: Literal["normal", "italic", "oblique"] = "normal",
    xticks_fontsize: float = 12.0,
    xticks_color: str = "black",
    xticks_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    xticks_style: Literal["normal", "italic", "oblique"] = "normal",
    yticks_fontsize: float = 12.0,
    yticks_color: str = "black",
    yticks_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    yticks_style: Literal["normal", "italic", "oblique"] = "normal",
    grid: bool = True,
    grid_linestyle: str = "--",
    grid_alpha: float = 0.7,
    background_color: str | None = "white",
    sort_order: Literal["ascending", "descending"] = "descending",
) -> KgnBarplot:
    """Generate a publication-grade customizable horizontal barplot for pathway completeness.

    Transforms eggNOG-mapper output columns into molten score matrices, drops non-numeric
    or empty reference pathways, normalizes categorical layout indexes, and maps a discrete
    gradient color palette across completion values.

    Args:
        df: Input DataFrame containing parsed pathway completeness vectors.
        figsize: Geometric allocation limits (width, height) defining canvas borders.
        cmap: Target string name lookup or a sequential list of direct hexadecimal colors.
        cmap_range: Truncation boundary offsets (start, end) utilized to sample the color map.
        title: Global text message identifier rendering above the drawing matrix.
        title_fontsize, title_color, title_weight, title_style: Font properties for title.
        xlabel, ylabel: Text content values mapping coordinates descriptors.
        xlabel_fontsize, ylabel_fontsize: Typography scale indices.
        xlabel_color, ylabel_color: Text color variables mapping target labels.
        xlabel_weight, ylabel_weight: Structural typographic density metrics.
        xlabel_style, ylabel_style: Geometric font slope configurations.
        xticks_fontsize, xticks_color, xticks_weight, xticks_style: X-tick typography rules.
        yticks_fontsize, yticks_color, yticks_weight, yticks_style: Y-tick typography rules.
        grid: Toggles background coordinate reference line structures.
        grid_linestyle: Grid line texture rendering parameter.
        grid_alpha: Opacity index managing visibility bounds of grid elements.
        background_color: Primary layout canvas backdrop color mapping.
        sort_order: Strict directional sort state applied across completeness metrics.

    Returns:
        KgnBarplot: Container instance holding references to optimized figures.
    """
    # Validate sample elements and apply strict ordered categorical indices
    working_df = df.drop(columns=["Function"], errors="ignore")
    working_df = working_df.select_dtypes(include="number")
    working_df = working_df.loc[:, (working_df > 0).any(axis=0)]
    df_melted = working_df.melt(var_name="Pathway", value_name="Score")

    df_melted = df_melted.sort_values(
        by="Score", ascending=(sort_order == "ascending")
    ).reset_index(drop=True)

    df_melted["Pathway"] = df_melted["Pathway"].apply(
        lambda x: x.capitalize() if isinstance(x, str) and x.islower() else x
    )

    # Establish palette map dictionaries compliant with static analysis
    n_bars = df_melted["Pathway"].nunique()

    if isinstance(cmap, str):
        palette = sns.color_palette(cmap, n_colors=cmap_range[1])[cmap_range[0] :]
        palette = palette[:n_bars]  # trim to actual bar count
    elif isinstance(cmap, (list, tuple)) or hasattr(cmap, "__iter__"):
        palette = list(cmap)[:n_bars]  # trim to actual bar count
    else:
        palette = sns.color_palette("Greens", n_colors=cmap_range[1])[cmap_range[0] :]
        palette = palette[:n_bars]  # trim to actual bar count

    # Initialize structural subplots container canvas layers
    fig, ax = plt.subplots(figsize=figsize, facecolor=background_color)

    sns.barplot(
        data=df_melted,
        x="Score",
        y="Pathway",
        hue="Score",
        palette=palette,
        dodge=False,
        legend=False,
        ax=ax,
    )

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

    # Configure ticks geometric parameters and close active canvas stream descriptors
    for label in ax.get_xticklabels():
        label.set_fontsize(xticks_fontsize)
        label.set_color(xticks_color)
        label.set_fontweight(xticks_weight)
        label.set_fontstyle(xticks_style)

    for label in ax.get_yticklabels():
        label.set_fontsize(yticks_fontsize)
        label.set_color(yticks_color)
        label.set_fontweight(yticks_weight)
        label.set_fontstyle(yticks_style)

    if grid:
        ax.grid(axis="x", linestyle=grid_linestyle, alpha=grid_alpha, zorder=0)

    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_axisbelow(True)

    plt.close(fig)

    return KgnBarplot(fig, ax)
