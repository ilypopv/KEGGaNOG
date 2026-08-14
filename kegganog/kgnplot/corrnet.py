"""Correlation network visualization module for KEGG module completion profiles.

This module builds standalone topological graph networks tracking sample-to-sample
reconstruction consistency and mathematical correlation strengths extracted from
eggNOG-mapper functional annotations.
"""

from typing import Literal

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import Colormap

from .base import KgnPlotBase


class KgnCorrnet(KgnPlotBase):
    """Orchestration context wrapper encapsulating Matplotlib correlation network layouts."""

    def __init__(self, fig: plt.Figure, ax: plt.Axes) -> None:
        """Initialize the correlation network canvas with layout metrics.

        Args:
            fig: The Matplotlib Figure container hosting the drawing canvas.
            ax: The core underlying Axes coordinate grid mapper.
        """
        super().__init__(fig, ax)


def correlation_network(
    df: pd.DataFrame,
    figsize: tuple[float, float] = (12.0, 6.0),
    threshold: float = 0.5,
    node_size: float = 700.0,
    node_color: str = "#A3D5FF",
    node_edgecolors: str = "#03045E",
    node_linewidths: float = 1.5,
    label_fontsize: float = 8.0,
    label_color: str = "#03045E",
    label_verticalalignment: Literal["center", "top", "bottom", "baseline"] = "center",
    label_horizontalalignment: Literal["center", "right", "left"] = "center",
    label_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    edge_cmap: str | Colormap = colormaps["coolwarm"],
    cbar_size: float = 0.5,
    title: str | None = None,
    title_fontsize: float = 16.0,
    title_color: str = "black",
    title_weight: Literal["normal", "bold", "heavy", "light"] = "normal",
    title_style: Literal["normal", "italic", "oblique"] = "normal",
    background_color: str | None = "white",
    save_matrix: str | None = None,
) -> KgnCorrnet:
    """Generate a publication-grade customizable correlation network for KEGG reconstructions.

    Computes cross-sample Pearson correlation matrices from functional profiles,
    filters topological edges using strict minimal thresholds, normalizes graph layout vectors,
    and maps geometric edge properties to statistical linkage strengths.

    Args:
        df: Input DataFrame containing parsed pathway completeness vectors.
        figsize: Geometric allocation limits (width, height) defining canvas borders.
        threshold: Absolute minimal correlation cut-off required to retain graph edge linkages.
        node_size: Normalized scalar volume index assigned to network nodes.
        node_color: Filling color identification token mapped to graph nodes.
        node_edgecolors: Border outline color identifier separating adjacent nodes.
        node_linewidths: Structural boundary border scale factor mapped onto nodes.
        label_fontsize: Typography size constraint allocated to sample text markers.
        label_color: Text color identifier allocated to central sample labels.
        label_verticalalignment: Geometric text layout vertical vector alignment constraint.
        label_horizontalalignment: Geometric text layout horizontal vector alignment constraint.
        label_weight: Structural typographic density metric specified for graph labels.
        edge_cmap: Target string colormap descriptor or a native Matplotlib Colormap object.
        cbar_size: Scaled metric index tracking colorbar widget dimensions.
        title: Global text message identifier rendering above the drawing matrix.
        title_fontsize, title_color, title_weight, title_style: Font properties for title.
        background_color: Primary layout canvas backdrop color mapping.
        save_matrix: Explicit destination file path to export computed cross-correlation TSV matrix.

    Returns:
        KgnCorrnet: Container instance holding references to optimized figures.
    """
    # Compute cross-correlation matrix from target numeric profiles
    correlation_matrix = df.iloc[:, 1:].corr()
    cor_threshold = threshold

    # Initialize topological graph object framework
    G = nx.Graph()

    for col in correlation_matrix.columns:
        G.add_node(col)

    # Populate network edges based on threshold filtration rules
    edges_list = []
    for i, col1 in enumerate(correlation_matrix.columns):
        for j, col2 in enumerate(correlation_matrix.columns):
            if i < j and abs(correlation_matrix.iloc[i, j]) > cor_threshold:
                weight = abs(correlation_matrix.iloc[i, j])
                edges_list.append((col1, col2, weight))

    G.add_weighted_edges_from(edges_list)

    # Calibrate dynamic edge width arrays to reflect linkage strength
    if G.number_of_edges() > 0:
        weights = [d["weight"] for _, _, d in G.edges(data=True)]
        max_weight = max(weights)
        min_weight = min(weights)

        edge_widths = [
            (
                (1.0 + (w - min_weight) / (max_weight - min_weight))
                if max_weight > min_weight
                else 2.0
            )
            for w in weights
        ]
    else:
        weights = []
        edge_widths = []

    # Initialize structural subplots container canvas layers
    fig, ax = plt.subplots(figsize=figsize, facecolor=background_color)
    pos = nx.spring_layout(G, seed=42)

    # Map geometric nodes, links, and label annotations onto target axis
    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=node_size,
        node_color=node_color,
        edgecolors=node_edgecolors,
        linewidths=node_linewidths,
        ax=ax,
    )

    drawn_edges = nx.draw_networkx_edges(
        G,
        pos,
        width=edge_widths,
        alpha=0.7,
        edge_color=weights,
        edge_cmap=edge_cmap,
        ax=ax,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        font_size=label_fontsize,
        font_color=label_color,
        verticalalignment=label_verticalalignment,
        horizontalalignment=label_horizontalalignment,
        font_weight=label_weight,
        ax=ax,
    )

    # Map colorbar indicators mapping weight vectors
    if G.number_of_edges() > 0 and drawn_edges is not None:
        cbar = plt.colorbar(drawn_edges, shrink=cbar_size, ax=ax)
        cbar.set_label("Correlation Strength")

    # Apply customized typography parameters to coordinate boundaries
    if title:
        ax.set_title(
            title,
            fontsize=title_fontsize,
            color=title_color,
            weight=title_weight,
            style=title_style,
        )

    ax.axis("off")

    # Export computed matrix workspace data if pathways are specified
    if save_matrix:
        correlation_matrix.to_csv(save_matrix, sep="\t")
        print(f"Correlation matrix saved as {save_matrix}")

    # Close active canvas stream descriptors to isolate thread states
    plt.close(fig)

    return KgnCorrnet(fig, ax)
