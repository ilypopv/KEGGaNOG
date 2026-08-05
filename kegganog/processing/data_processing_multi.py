"""Data parsing, execution, and matrix aggregation engine for multi-cohort workflows.

This module orchestrates parallel or sequential functional evaluation of multiple
eggNOG-mapper runs, handles matrix transposition, and merges individual pathway
reconstructions into unified cohort comparative tables.
"""

from __future__ import annotations

import csv
import io
import logging
import subprocess
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from .data_processing import _ensure_kegg_decoder

# Initialize module-level isolated logger
_log: logging.Logger = logging.getLogger(__name__)


def parse_emapper(
    input_file: Path | str, sample_folder: Path | str, file_prefix: str
) -> Path:
    """Extract and format a single eggNOG annotation table within a multi-sample cohort run.

    Args:
        input_file: Target raw input path containing annotation tables.
        sample_folder: Dedicated sample subdirectory allocated for temporary buffers.
        file_prefix: Unique sample identifier string used for isolation labeling.

    Returns:
        Path: Location pointer to the reformatted cohort-ready KO map table.

    Raises:
        KeyError: If the mandatory structural 'KEGG_ko' metadata column is absent.
    """
    in_path: Path = Path(input_file)
    sample_path: Path = Path(sample_folder)

    _log.info("Parsing eggNOG-mapper source file: %s", in_path)

    # Scan annotation document lines securely to detect active header tables
    with open(in_path, encoding="utf-8") as fh:
        header_row: int = next(
            i for i, line in enumerate(fh) if line.strip() and not line.startswith("##")
        )

    # Read target segments into memory using Pandas with live progress metrics
    with tqdm(total=1, desc=f"Reading {file_prefix}") as pbar:
        df_filtered: pd.DataFrame = pd.read_csv(
            in_path, sep="\t", skiprows=header_row, header=0
        )
        pbar.update(1)

    # Enforce operational schema integrity check
    if "KEGG_ko" not in df_filtered.columns:
        raise KeyError(
            f"Mandatory 'KEGG_ko' column not found in {in_path.name}. Please verify file format integrity."
        )

    # Extract valid assigned functional metrics block
    df_kegg_ko: pd.DataFrame = df_filtered[["KEGG_ko"]].copy()
    df_kegg_ko = df_kegg_ko[df_kegg_ko["KEGG_ko"] != "-"]

    # Format specific KO map column spaces using the isolated cohort prefix keys
    with tqdm(total=1, desc=f"Formatting {file_prefix}") as pbar:
        df_kegg_ko["KEGG_ko"] = df_kegg_ko["KEGG_ko"].str.replace(
            r"ko:(K\d+)", rf"{file_prefix} \1", regex=True
        )
        df_kegg_ko["KEGG_ko"] = df_kegg_ko["KEGG_ko"].str.replace(",", "\n")

        # Pipe string buffers inside system memory using IO representations
        buffer: io.StringIO = io.StringIO()
        df_kegg_ko.to_csv(
            buffer,
            sep="\t",
            index=False,
            header=False,
            quoting=csv.QUOTE_MINIMAL,
            escapechar="\\",
        )
        buffer.seek(0)

        content: str = buffer.read().replace('"', "")
        pbar.update(1)

    # Flush text stream metrics out to physical disk spaces
    parsed_filtered_file: Path = sample_path / f"{file_prefix}_parsed_KO_terms.txt"
    parsed_filtered_file.write_text(content, encoding="utf-8")

    return parsed_filtered_file


def run_kegg_decoder(
    input_file: Path | str, sample_folder: Path | str, file_prefix: str
) -> Path:
    """Execute KEGG-Decoder runner for a specific cohort sample matrix slice.

    Args:
        input_file: Source path pointing to the pre-parsed KO map documents.
        sample_folder: Target location container directory for this specific sample.
        file_prefix: Unique target sample identity prefix.

    Returns:
        Path: Location pointer to the generated output pathway TSV table.
    """
    in_path: Path = Path(input_file)
    sample_path: Path = Path(sample_folder)
    output_file: Path = sample_path / f"{file_prefix}_pathways.tsv"

    _log.info("Executing KEGG-Decoder sub-process wrapper on sample: %s", file_prefix)

    # Locate and validate background engine script dependencies
    package_dir: Path = Path(__file__).resolve().parent
    kegg_decoder_script: Path = package_dir / "KEGG_decoder.py"
    _ensure_kegg_decoder(kegg_decoder_script)

    # Trigger subprocess runtime pipeline call
    with tqdm(total=1, desc="Executing KEGG-Decoder") as pbar:
        command: list[str] = [
            "python",
            str(kegg_decoder_script),
            "-i",
            str(in_path),
            "-o",
            str(output_file),
        ]
        subprocess.run(command, check=True)
        pbar.update(1)

    # Perform text transformations to fix column metadata layouts
    lines: list[str] = output_file.read_text(encoding="utf-8").splitlines(keepends=True)

    if lines:
        first_line: list[str] = lines[0].strip().split("\t")
        processed_first_line: list[str] = [
            x.capitalize() if isinstance(x, str) and x.islower() else x
            for x in first_line
        ]
        lines[0] = "\t".join(processed_first_line) + "\n"

    # Evict structural placeholder labels in favor of active sample prefix tags
    content: str = "".join(lines).replace("SAMPLE", file_prefix)
    output_file.write_text(content, encoding="utf-8")

    return output_file


def merge_outputs(output_folder: Path | str) -> pd.DataFrame:
    """Aggregate individual cohort pathway matrix slices into a single master TSV document.

    Transposes matrix tracks, matches functional indexes, and outputs an outer-joined
    comparative table across all discovered valid sample blocks.

    Args:
        output_folder: Root output directory hosting the target 'temp_files' space.

    Returns:
        pd.DataFrame: Merged cohort analytical data matrix layout.
    """
    out_dir: Path = Path(output_folder)
    _log.info("Aggregating individual multi-sample pathway profiles...")

    merged_df: pd.DataFrame = pd.DataFrame()

    # Scan target directories securely using clean Path.glob expressions
    temp_dir: Path = out_dir / "temp_files"
    output_files: list[Path] = list(temp_dir.glob("*/*_pathways.tsv"))

    # Iterate and parse metrics blocks into transposable structures
    for file_path in output_files:
        file_prefix: str = file_path.stem.replace("_pathways", "")
        df: pd.DataFrame = pd.read_csv(file_path, sep="\t", index_col=0)

        df_transposed: pd.DataFrame = df.T
        df_transposed.columns = [file_prefix]

        if merged_df.empty:
            merged_df["Function"] = df_transposed.index

        merged_df = pd.merge(
            merged_df, df_transposed, left_on="Function", right_index=True, how="outer"
        )

    # Flush final integrated dataset to a structural TSV on disk
    merged_output_file: Path = out_dir / "merged_pathways.tsv"
    merged_df.to_csv(merged_output_file, sep="\t", index=False)
    _log.info(
        "Cohort files aggregated successfully into master matrix: '%s'",
        merged_output_file,
    )

    return merged_df
