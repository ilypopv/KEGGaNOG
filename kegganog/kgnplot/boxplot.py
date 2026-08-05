"""Boxplot visualization module for KEGG module completion profiles.

This module builds standalone box-and-whisker plots tracking macro-distribution
metrics and pathway completeness density alterations across independent samples.
"""

from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .base import KgnPlotBase


class KgnBoxplot(KgnPlotBase):
    """Orchestration context wrapper encapsulating Matplotlib boxplot layouts."""

    def __init__(self, fig: plt.Figure, ax: plt.Axes) -> None:
        """Initialize the boxplot canvas with layout metrics.

        Args:
            fig: The Matplotlib Figure container hosting the drawing canvas.
            ax: The core underlying Axes coordinate grid mapper.
        """
        super().__init__(fig, ax)


def boxplot(
    df: pd.DataFrame,
    figsize: tuple[float, float] = (12.0, 6.0),
    color: str | None = "blue",
    showfliers: bool = True,
    title: str | None = None,
    title_fontsize: float = 16.0,
    title_color: str = "black",
    title_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    title_style: Literal["normal", "italic", "oblique"] = "normal",
    xlabel: str = "Samples",
    xlabel_fontsize: float = 14.0,
    xlabel_color: str = "black",
    xlabel_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    xlabel_style: Literal["normal", "italic", "oblique"] = "normal",
    ylabel: str = "Completeness Value",
    ylabel_fontsize: float = 14.0,
    ylabel_color: str = "black",
    ylabel_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    ylabel_style: Literal["normal", "italic", "oblique"] = "normal",
    xticks_rotation: float = 45.0,
    xticks_ha: Literal["left", "right", "center"] = "center",
    xticks_fontsize: float = 12.0,
    xticks_color: str = "black",
    xticks_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    xticks_style: Literal["normal", "italic", "oblique"] = "normal",
    yticks_fontsize: float = 12.0,
    yticks_color: str = "black",
    yticks_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    yticks_style: Literal["normal", "italic", "oblique"] = "normal",
    grid: bool = True,
    grid_color: str = "gray",
    grid_linestyle: str = "--",
    grid_linewidth: float = 0.5,
    background_color: str | None = "white",
) -> KgnBoxplot:
    """Generate a publication-grade customizable boxplot for pathway completeness distributions.

    Parses multi-sample reconstruction matrices, extracts numeric functional metrics,
    isolates coordinate matrices, and renders statistical quartile boxes profiling variation.

    Args:
        df: Input DataFrame containing parsed pathway completeness vectors.
        figsize: Geometric allocation limits (width, height) defining canvas borders.
        color: Primary uniform filling color hex identifier mapped to individual boxes.
        showfliers: Flag indicating whether statistical outlier points should be plotted.
        title: Global text message identifier rendering above the drawing matrix.
        title_fontsize, title_color, title_weight, title_style: Font properties for title.
        xlabel, ylabel: Text content values mapping coordinates descriptors.
        xlabel_fontsize, ylabel_fontsize: Typography scale indices.
        xlabel_color, ylabel_color: Text color variables mapping target labels.
        xlabel_weight, ylabel_weight: Structural typographic density metrics.
        xlabel_style, ylabel_style: Geometric font slope configurations.
        xticks_rotation, xticks_ha: Position variables mapping target X tick attributes.
        xticks_fontsize, xticks_color, xticks_weight, xticks_style: X-tick typography rules.
        yticks_fontsize, yticks_color, yticks_weight, yticks_style: Y-tick typography rules.
        grid: Toggles background coordinate reference line structures.
        grid_color: Target color identification hex parameter allocated to grid lines.
        grid_linestyle: Grid line texture rendering parameter.
        grid_linewidth: Linear density configuration of background alignment segments.
        background_color: Primary layout canvas backdrop color mapping.

    Returns:
        KgnBoxplot: Container instance holding references to optimized figures.
    """
    # Initialize structural subplots container canvas layers
    fig, ax = plt.subplots(figsize=figsize, facecolor=background_color)

    # Map underlying statistical matrices bypassing function label features
    sns.boxplot(data=df.iloc[:, 1:], color=color, showfliers=showfliers, ax=ax)

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
    plt.xticks(
        rotation=xticks_rotation,
        ha=xticks_ha,
        fontsize=xticks_fontsize,
        color=xticks_color,
        weight=xticks_weight,
        style=xticks_style,
    )

    plt.yticks(
        fontsize=yticks_fontsize,
        color=yticks_color,
        weight=yticks_weight,
        style=yticks_style,
    )

    if grid:
        ax.grid(color=grid_color, linestyle=grid_linestyle, linewidth=grid_linewidth)

    plt.close(fig)
    return KgnBoxplot(fig, ax)
