"""Radar plot visualization module for KEGG module completion profiles.

This module builds standalone polar spider/radar charts tracking individual
pathway completeness metrics across continuous sample coordinates.
"""

import warnings
from collections.abc import Sequence
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .base import KgnPlotBase


class KgnRadar(KgnPlotBase):
    """Orchestration context wrapper encapsulating Matplotlib polar radar layouts."""

    def __init__(self, fig: plt.Figure, ax: plt.Axes) -> None:
        """Initialize the radar plot canvas with layout metrics.

        Args:
            fig: The Matplotlib Figure container hosting the drawing canvas.
            ax: The core underlying Axes coordinate grid mapper.
        """
        super().__init__(fig, ax)


def radarplot(
    df: pd.DataFrame,
    pathways: Sequence[str],
    figsize: tuple[float, float] = (8.0, 8.0),
    colors: Sequence[str] | None = None,
    sample_order: Sequence[str] | None = None,
    title: str | None = None,
    title_fontsize: float = 14.0,
    title_color: str = "black",
    title_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    title_style: Literal["normal", "italic", "oblique"] = "normal",
    title_y: float = 1.1,
    label_fontsize: float = 10.0,
    label_color: str = "black",
    label_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    label_style: Literal["normal", "italic", "oblique"] = "normal",
    label_background: str | None = None,
    label_edgecolor: str | None = None,
    label_pad: float = 1.05,
    ytick_fontsize: float = 8.0,
    ytick_color: str = "black",
    ytick_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    ytick_alpha: float = 0.5,
    yticklabels: Sequence[str] | None = None,
    fill_alpha: float = 0.25,
    line_width: float = 2.0,
    line_style: Literal["solid", "dashed", "dashdot", "dotted", "-"] = "solid",
    background_color: str | None = "white",
    legend_loc: str = "upper right",
    legend_bbox: tuple[float, float] = (1.3, 1.1),
    show_legend: bool = True,
) -> KgnRadar:
    """Generate a publication-grade customizable polar radar chart for KEGG pathways.

    Extracts completeness vector rows matching target functional descriptions,
    maps coordinate indices to a radial layout geometry, project data shapes,
    and applies custom transparency blending keys across multiple targets.

    Args:
        df: Input DataFrame containing parsed pathway completeness vectors.
        pathways: Target collection containing pathway identifier strings (max 4).
        figsize: Geometric allocation limits (width, height) defining canvas borders.
        colors: Target sequential list containing custom hex/string colors for plotting.
        sample_order: Explicit layout sequence locking the display order of sample axes.
        title: Global text message identifier rendering above the drawing matrix.
        title_fontsize, title_color, title_weight, title_style: Font properties for title.
        title_y: Scaled vertical offset position parameter handling title layouts.
        label_fontsize, label_color, label_weight, label_style: Typography configuration.
        label_background: Background color identifier used to coat bounding label cells.
        label_edgecolor: Outline boundary color constraint assigned to textual label boxes.
        label_pad: Geometric distance buffer between radial boundaries and sample labels.
        ytick_fontsize, ytick_color, ytick_weight, ytick_alpha: Y-tick typography rules.
        yticklabels: Explicit custom list of labels applied directly to concentric grid paths.
        fill_alpha: Opacity scale constraint applied to drawing layers.
        line_width, line_style: Structural line property parameters mapped onto shapes.
        background_color: Primary layout canvas backdrop color mapping.
        legend_loc: Positional anchoring code identifier tracking layout widgets.
        legend_bbox: Coordinate anchor box offsets defining bounding regions for legends.
        show_legend: If False, completely suppresses widget layer execution.

    Returns:
        KgnRadar: Container instance holding references to optimized figures.

    Raises:
        ValueError: Triggered if requested pathway list array exceeds the limit of 4 items.
    """
    # Enforce mathematical constraints bounding active pathways array
    if len(pathways) > 4:
        raise ValueError("Maximum of 4 pathways can be plotted at once.")

    # Validate sample elements and apply strict ordered categorical indices
    if sample_order is None:
        sample_order = [col for col in df.columns if col != "Function"]

    # Reshape layout spreadsheet structure via pivot operations
    num_vars = len(sample_order)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]

    # Initialize structural subplots container canvas layers
    fig, ax = plt.subplots(
        figsize=figsize, subplot_kw={"polar": True}, facecolor=background_color
    )

    # Establish palette map dictionaries compliant with static analysis
    if colors is None:
        colors_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    else:
        colors_list = list(colors)

    # Map underlying statistical matrices bypassing function label features
    for i, function in enumerate(pathways):
        row = df[df["Function"] == function]
        if row.empty:
            warnings.warn(
                f"Pathway '{function}' not found in DataFrame. Skipping.",
                UserWarning,
                stacklevel=2,
            )
            continue
        values = row[sample_order].values.flatten().tolist()
        values += values[:1]

        ax.plot(
            angles,
            values,
            label=function,
            linewidth=line_width,
            linestyle=line_style,
            color=colors_list[i % len(colors_list)],
        )
        ax.fill(
            angles,
            values,
            color=colors_list[i % len(colors_list)],
            alpha=fill_alpha,
        )

    ax.set_xticks([])

    # Manually extrapolate background spokes mapping coordinate grid lines
    for angle in angles[:-1]:
        ax.plot(
            [angle, angle],
            [0, 1],
            color="lightgray",
            linewidth=1,
            linestyle="solid",
            zorder=0,
        )

    for angle, label in zip(angles[:-1], sample_order):
        ax.text(
            angle,
            label_pad,
            label,
            ha="center",
            va="center",
            fontsize=label_fontsize,
            color=label_color,
            fontweight=label_weight,
            style=label_style,
            bbox=(
                {
                    "facecolor": label_background if label_background else "none",
                    "edgecolor": label_edgecolor if label_edgecolor else "none",
                    "boxstyle": "round,pad=0.2",
                }
                if label_background or label_edgecolor
                else None
            ),
        )

    # Configure ticks geometric parameters and close active canvas stream descriptors
    grid_vals = np.linspace(0.2, 1.0, 5)
    ax.set_yticks(grid_vals)
    ax.set_ylim(0, 1.0)

    if yticklabels is None:
        calculated_labels = [""] * len(grid_vals)
        calculated_labels[grid_vals.tolist().index(0.2)] = "0.2"
        calculated_labels[grid_vals.tolist().index(1.0)] = "1.0"
    else:
        calculated_labels = list(yticklabels)

    ax.set_yticklabels(
        calculated_labels,
        fontsize=ytick_fontsize,
        color=ytick_color,
        fontweight=ytick_weight,
        alpha=ytick_alpha,
    )

    # Apply customized typography parameters to coordinate boundaries
    if title:
        plt.title(
            title,
            size=title_fontsize,
            color=title_color,
            weight=title_weight,
            style=title_style,
            y=title_y,
        )

    # Build legend overlays if requested by execution flags
    if show_legend:
        ax.legend(loc=legend_loc, bbox_to_anchor=legend_bbox)

    plt.close(fig)

    return KgnRadar(fig, ax)
