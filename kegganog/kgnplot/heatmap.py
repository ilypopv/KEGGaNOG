"""Universal heatmap orchestration module for KEGG module completion profiles.

This module acts as a unified factory interface routing data profiles between
single-sample and multi-sample matrix heatmaps, managing runtime isolation
and file system cleanup.
"""

import contextlib
import io
import shutil
import tempfile
import warnings
from collections.abc import Generator, Sequence
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import tqdm

# Fix: Use explicit relative package imports to guarantee internal module visibility
from ..cheatmaps.grouped_heatmap import generate_grouped_heatmap
from ..cheatmaps.grouped_heatmap_multi import generate_grouped_heatmap_multi
from ..cheatmaps.simple_heatmap import generate_heatmap
from ..cheatmaps.simple_heatmap_multi import generate_heatmap_multi
from .base import KgnPlotBase


@contextlib.contextmanager
def silent_plot_and_tqdm() -> Generator[None, None, None]:
    """Context manager to completely suppress console output, tqdm bars, and UI renders."""

    original_show = plt.show
    original_tqdm = tqdm.tqdm

    def silent_show(*args: Any, **kwargs: Any) -> None:
        pass

    plt.show = silent_show

    class SilentTqdm(tqdm.tqdm):
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            kwargs["disable"] = True
            super().__init__(*args, **kwargs)

    tqdm.tqdm = SilentTqdm

    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        try:
            yield
        finally:
            plt.show = original_show
            tqdm.tqdm = original_tqdm


class KgnHeatmap(KgnPlotBase):
    """Orchestration context wrapper encapsulating Matplotlib heatmap canvas layouts."""

    def __init__(self, fig: plt.Figure, ax: Sequence[plt.Axes]) -> None:
        """Initialize the heatmap canvas with layout metrics.

        Args:
            fig: The Matplotlib Figure container hosting the drawing canvas.
            ax: The core underlying Axes coordinate grid mapper.
        """
        super().__init__(fig, ax[0])
        self.ax: Sequence[plt.Axes] = ax


def heatmap(
    df: pd.DataFrame,
    figsize: tuple[float, float] | None = None,
    color: str | None = None,
    group: bool = False,
    sample_name: str | None = None,
    annot: bool = True,
) -> KgnHeatmap:
    """Universal heatmap router factory for single and multi-sample KEGG reconstructions.

    Serializes runtime DataFrames into structural files, routes matrices to dedicated
    sub-package plotting targets, dynamically configures categorical groupings,
    and wraps output figure assets securely.

    Args:
        df: Input DataFrame containing parsed pathway completeness vectors.
        figsize: Geometric allocation limits (width, height) defining canvas borders.
        color: Target string colormap descriptor (e.g., "Greens", "Blues").
        group: Flag indicating whether layout modules should group pathways by category.
        sample_name: Optional sample name label override used exclusively in single-mode.
        annot: Toggles text value annotations inside heatmaps cells (ignored in multi-mode).

    Returns:
        KgnHeatmap: Container instance holding references to optimized figures.

    Raises:
        ValueError: Triggered if the input matrix contains empty rows or invalid geometry.
    """
    # Serialize runtime matrix workspace files into secure temporary pathways
    with tempfile.NamedTemporaryFile(delete=False, mode="w", newline="") as temp_file:
        df.to_csv(temp_file, sep="\t", index=False)
        temp_file_path = temp_file.name

    # Establish default fallback color mapping profiles
    target_color = color or "Blues"

    is_single = df.shape[0] == 1
    is_multi = df.shape[0] > 1

    if not (is_single or is_multi):
        raise ValueError(
            "DataFrame does not fit the expected dimensions for heatmap generation"
        )

    # Resolve and execute target rendering engine using explicit branch isolation
    output_folder = tempfile.mkdtemp()

    with silent_plot_and_tqdm():
        if not is_single:
            if sample_name is not None or not annot:
                warnings.warn(
                    "The 'sample_name' and 'annot' arguments are ignored for multi-heatmaps.",
                    UserWarning,
                    stacklevel=2,
                )

            if group:
                fig, ax = generate_grouped_heatmap_multi(
                    kegg_decoder_file=df,
                    output_folder=output_folder,
                    dpi=300,
                    color=target_color,
                    figsize=figsize,
                )
            else:
                fig, ax = generate_heatmap_multi(
                    kegg_decoder_file=df,
                    output_folder=output_folder,
                    dpi=300,
                    color=target_color,
                    figsize=figsize,
                )
        else:
            # Single-sample mode branches
            if group:
                fig, ax = generate_grouped_heatmap(
                    kegg_decoder_file=temp_file_path,
                    output_folder=output_folder,
                    dpi=300,
                    color=target_color,
                    sample_name=sample_name,
                    figsize=figsize,
                    annot=annot,
                )
            else:
                fig, ax = generate_heatmap(
                    kegg_decoder_file=temp_file_path,
                    output_folder=output_folder,
                    dpi=300,
                    color=target_color,
                    sample_name=sample_name,
                    figsize=figsize,
                    annot=annot,
                )

    # Perform safe structural system purge cleaning temporary layout buffers
    shutil.rmtree(output_folder, ignore_errors=True)

    # Close canvas stream descriptors to completely isolate thread states
    plt.close(fig)

    return KgnHeatmap(fig, ax)
