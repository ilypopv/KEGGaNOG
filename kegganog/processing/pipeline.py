"""Central orchestration hub managing coordinate execution tracks for KEGGaNOG.

This module exposes unified interface wrappers to execute single-sample or
multi-sample comparative data extraction flows, handles transient memory files,
routes tracking parameters, and bundles cross-platform layout packages.
"""

from __future__ import annotations

import contextlib
import logging
import os
import shutil
import tempfile
import zipfile
from collections.abc import Generator
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ..schemas import ValidColor
from .data_processing import parse_emapper, run_kegg_decoder
from .data_processing_multi import merge_outputs
from .data_processing_multi import parse_emapper as parse_emapper_multi
from .data_processing_multi import run_kegg_decoder as run_kegg_decoder_multi

# Initialize module-level isolated logger
_log: logging.Logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Structural Runtime Containers
# ---------------------------------------------------------------------------


@dataclass
class PipelineResult:
    """Core container hosting location parameters for runtime pipeline assets.

    Attributes:
        zip_path: Stable path pointer targeting packaged result archives.
        png_path: Target layout pointer referencing the primary matrix visualization.
        tsv_path: Location pointer targeting the structured summary matrix data sheet.
        samples: Array sequence containing verified sample identity keys.
        pathways: Sequence hosting discovered functional pathway identifiers.
    """

    zip_path: Path
    png_path: Path
    tsv_path: Path
    samples: list[str] = field(default_factory=list)
    pathways: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Internal Helpers & Packed Asset Handlers
# ---------------------------------------------------------------------------


def _secure_temp_path(suffix: str) -> Path:
    """Acquire a securely closed, atomic temporary system file pointer.

    Args:
        suffix: File extension descriptor type assignment.

    Returns:
        Path: Valid closed location reference pointer inside OS temp storage.
    """
    fd, path = tempfile.mkstemp(suffix=suffix)

    os.close(fd)
    return Path(path)


def _pack_results(output_dir: Path) -> tuple[Path, Path]:
    """Compress structural result files into centralized deployment package files.

    Args:
        output_dir: Base directory zone targeting completed profile tables.

    Returns:
        tuple[Path, Path]: Paired path tokens for the generated (.zip, .png) assets.

    Raises:
        FileNotFoundError: If evaluation sequences fail to generate matrix charts.
    """
    png_files: list[Path] = list(output_dir.rglob("*.png"))

    if not png_files:
        raise FileNotFoundError(
            "Analytical run finished successfully but failed to output visual chart matrices. "
            "Please verify that your input files contain valid orthology-based functional maps."
        )

    stable_png: Path = _secure_temp_path(".png")
    shutil.copy2(png_files[0], stable_png)

    stable_zip: Path = _secure_temp_path(".zip")
    with zipfile.ZipFile(stable_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for result_file in output_dir.rglob("*"):
            if result_file.is_file():
                zf.write(result_file, result_file.relative_to(output_dir))

    return stable_zip, stable_png


@contextlib.contextmanager
def _resolve_output_context(
    output_dir: Path | str | None,
) -> Generator[Path, None, None]:
    """Yield a stable active target directory context matching CLI or Web runstates.

    Args:
        output_dir: User-specified destination path pointer if in CLI execution mode.

    Yields:
        Path: Valid workspace folder path context.
    """
    if output_dir is not None:
        out_path: Path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        yield out_path
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "output"
            out_path.mkdir()
            yield out_path


# ---------------------------------------------------------------------------
# Core Background Execution Logic Layers
# ---------------------------------------------------------------------------


def _run_single_in_dir(
    file_bytes: bytes,
    sample_name: str,
    dpi: int,
    color: ValidColor,
    group: bool,
    output_dir: Path,
) -> tuple[list[str], list[str], Path]:
    """Execute single-sample functional table mapping extraction matrix layouts."""
    from ..cheatmaps import grouped_heatmap, simple_heatmap

    temp_folder: Path = output_dir / "temp_files"
    temp_folder.mkdir(exist_ok=True)

    input_file: Path = temp_folder / "input.annotations"
    input_file.write_bytes(file_bytes)

    # Decode target sequences utilizing abstracted translation engine components
    profile_file: Path = run_kegg_decoder(
        parse_emapper(input_file, temp_folder),
        output_dir,
        sample_name,
    )

    if group:
        grouped_heatmap.generate_grouped_heatmap(
            str(profile_file), str(output_dir), dpi, color, sample_name
        )
    else:
        simple_heatmap.generate_heatmap(
            str(profile_file), str(output_dir), dpi, color, sample_name
        )

    df: pd.DataFrame = pd.read_csv(profile_file, sep="\t", index_col=0)
    pathways: list[str] = list(df.columns)

    return [sample_name], pathways, profile_file


def _run_multi_in_dir(
    named_files: list[tuple[str, bytes]],
    dpi: int,
    color: ValidColor,
    group: bool,
    output_dir: Path,
) -> tuple[list[str], list[str], Path]:
    """Execute multi-sample comparative cohort evaluation matrices routines."""
    from ..cheatmaps import grouped_heatmap_multi, simple_heatmap_multi

    temp_folder: Path = output_dir / "temp_files"
    temp_folder.mkdir(exist_ok=True)

    for filename, file_bytes in named_files:
        file_prefix: str = filename.replace(".emapper.annotations", "")
        if "." in file_prefix:
            file_prefix = Path(file_prefix).stem

        sample_folder: Path = temp_folder / file_prefix
        sample_folder.mkdir(exist_ok=True)

        input_file: Path = sample_folder / filename
        input_file.write_bytes(file_bytes)

        parsed_file: Path = parse_emapper_multi(input_file, sample_folder, file_prefix)
        run_kegg_decoder_multi(parsed_file, sample_folder, file_prefix)

    merged_df: pd.DataFrame = merge_outputs(output_dir)

    if group:
        grouped_heatmap_multi.generate_grouped_heatmap_multi(
            merged_df, str(output_dir), dpi, color
        )
    else:
        simple_heatmap_multi.generate_heatmap_multi(
            merged_df, str(output_dir), dpi, color
        )

    samples: list[str] = [c for c in merged_df.columns if c != "Function"]
    pathways: list[str] = list(merged_df["Function"].dropna().unique())
    merged_tsv: Path = output_dir / "merged_pathways.tsv"

    return samples, pathways, merged_tsv


# ---------------------------------------------------------------------------
# Public Unified Execution Gateways
# ---------------------------------------------------------------------------


def run_single(
    file_bytes: bytes,
    sample_name: str,
    dpi: int,
    color: ValidColor,
    group: bool,
    output_dir: Path | str | None = None,
) -> PipelineResult:
    """Execute the full core functional annotation profile pathway pipeline.

    Args:
        file_bytes: Raw text byte content arrays read from eggNOG source tables.
        sample_name: Explicit string tag assigned for column index labels.
        dpi: Pixel configuration layout multiplier applied to matrix charts.
        color: Target valid seaborn color map palette identifier literal.
        group: Flag forcing hierarchical visual grouping rows layout.
        output_dir: Explicit folder directory pointer (CLI mode token), or None
            to instruct dynamic temporal allocation (Web server UI context mode).

    Returns:
        PipelineResult: Standard data object referencing assets location markers.
    """
    with _resolve_output_context(output_dir) as out:
        samples, pathways, profile_file = _run_single_in_dir(
            file_bytes, sample_name, dpi, color, group, out
        )
        stable_tsv: Path = _secure_temp_path(".tsv")
        shutil.copy2(profile_file, stable_tsv)
        zip_path, png_path = _pack_results(out)

    return PipelineResult(
        zip_path=zip_path,
        png_path=png_path,
        tsv_path=stable_tsv,
        samples=samples,
        pathways=pathways,
    )


def run_multi(
    named_files: list[tuple[str, bytes]],
    dpi: int,
    color: ValidColor,
    group: bool,
    output_dir: Path | str | None = None,
) -> PipelineResult:
    """Execute the aggregate multi-sample cohort analytical matrix workflow.

    Args:
        named_files: Sequence hosting pairs of (filename_string, content_bytes).
        dpi: Target evaluation layout plot graphic density index.
        color: Valid Seaborn palette mapping assignment token.
        group: Toggle routing conditional categorization row structures.
        output_dir: Target destination workspace root path context pointer.

    Returns:
        PipelineResult: Consistent metadata structure hosting exported paths tokens.
    """
    with _resolve_output_context(output_dir) as out:
        samples, pathways, merged_tsv = _run_multi_in_dir(
            named_files, dpi, color, group, out
        )
        stable_tsv: Path = _secure_temp_path(".tsv")
        shutil.copy2(merged_tsv, stable_tsv)
        zip_path, png_path = _pack_results(out)

    return PipelineResult(
        zip_path=zip_path,
        png_path=png_path,
        tsv_path=stable_tsv,
        samples=samples,
        pathways=pathways,
    )
