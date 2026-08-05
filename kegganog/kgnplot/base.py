"""Base lifecycle management module for KEGGaNOG mathematical graphics.

This module establishes the core architectural contract for handling matplotlib
Figure and Axes lifecycles, preventing memory leaks in automated pipelines
by wrapping heavy graphical contexts into explicit transaction objects.
"""

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt


class KgnPlotBase:
    """Base transaction wrapper for Matplotlib figure lifecycle and resource isolation."""

    def __init__(self, fig: plt.Figure, ax: plt.Axes) -> None:
        """Initialize the graphical context container.

        Args:
            fig: Matplotlib Figure instance to handle.
            ax: Matplotlib Axes instance associated with the figure.
        """
        self.fig: plt.Figure = fig
        self.ax: plt.Axes = ax

    def plotfig(self) -> plt.Figure:
        """Finalize figure layout boundaries without invoking blocking GUI loops.

        Returns:
            plt.Figure: The standalone finalized figure object decoupled from internal lifecycle.
        """
        self.fig.tight_layout()
        return self.fig

    def savefig(
        self,
        path: str | Path,
        dpi: int = 300,
        transparent: bool = False,
        bbox_inches: Literal["tight"] | str | None = "tight",
    ) -> None:
        """Save the enclosed figure to a file using cross-platform safe IO pathways.

        Args:
            path: Destination file path (supports string or cross-platform pathlib.Path objects).
            dpi: Dots per inch resolution of the exported graphic matrix. Defaults to 300.
            transparent: Flag to determine background alpha transparency state. Defaults to False.
            bbox_inches: Bounding box configuration or literal "tight" for strict edge clipping.
                Defaults to "tight".
        """
        # Convert path argument to standard pathlib object for cross-platform robustness
        target_path: Path = Path(path)

        # Delegate byte transaction to Matplotlib engine with explicit settings
        self.fig.savefig(
            target_path,
            dpi=dpi,
            transparent=transparent,
            bbox_inches=bbox_inches,
        )
