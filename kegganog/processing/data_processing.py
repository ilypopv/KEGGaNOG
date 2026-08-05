"""Data parsing and execution sub-engine for internal functional profile workflows.

This module handles the extraction and formatting of KO (KEGG Orthology) terms
from eggNOG-mapper annotation source tables, manages dynamic on-demand runtime
downloads of the core translation script, and executes subprocess wrappers.
"""

from __future__ import annotations

import csv
import io
import logging
import subprocess
from pathlib import Path

import pandas as pd
from tqdm import tqdm

# Initialize module-level isolated logger
_log: logging.Logger = logging.getLogger(__name__)

_SCRIPT_URL: str = "https://raw.githubusercontent.com/bjtully/BioData/master/KEGGDecoder/KEGG_decoder.py"
_REQUEST_TIMEOUT: int = 10


def _ensure_kegg_decoder(script_path: Path) -> None:
    """Download the third-party KEGG_decoder.py dependencies on first use.

    Args:
        script_path: Targeted destination file path container.

    Raises:
        RuntimeError: Wrapper fault message if network IO operations fail.
    """
    if script_path.exists():
        return

    _log.info(
        "KEGG_decoder.py script dependency absent. Automating download from:\n  %s",
        _SCRIPT_URL,
    )
    try:
        import requests

        response: requests.Response = requests.get(
            _SCRIPT_URL, timeout=_REQUEST_TIMEOUT
        )
        response.raise_for_status()
        script_path.write_bytes(response.content)
        _log.info(
            "Successfully deployed KEGG_decoder.py dependency (%d bytes).",
            script_path.stat().st_size,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Failed to download KEGG_decoder.py: {exc}\n"
            f"Please download it manually from:\n  {_SCRIPT_URL}\n"
            f"and place it directly at target destination:  {script_path}"
        ) from exc


def parse_emapper(input_file: Path | str, temp_folder: Path | str) -> Path:
    """Extract and reformat raw eggNOG functional annotations into explicit KO maps.

    Args:
        input_file: Target raw input path containing annotation tables.
        temp_folder: Temporary runtime folder allocated for storage buffers.

    Returns:
        Path: Target location pointer to the newly reformatted KO map table.
    """
    in_path: Path = Path(input_file)
    temp_path: Path = Path(temp_folder)

    # Scan annotation document lines to detect active header tables
    with open(in_path, encoding="utf-8") as fh:
        header_row: int = next(
            i for i, line in enumerate(fh) if line.strip() and not line.startswith("##")
        )

    # Read data segments into Pandas DataFrame with live tracker metrics
    with tqdm(total=1, desc="Reading eggNOG-mapper annotations") as pbar:
        df_filtered: pd.DataFrame = pd.read_csv(
            in_path, sep="\t", skiprows=header_row, header=0
        )
        pbar.update(1)

    # Extract and drop missing or unassigned functional metrics
    df_kegg_ko: pd.DataFrame = df_filtered[["KEGG_ko"]].copy()
    df_kegg_ko = df_kegg_ko[df_kegg_ko["KEGG_ko"] != "-"]

    # Restructure annotation matrices layout into KEGG-Decoder target patterns
    with tqdm(total=2, desc="Formatting KEGG_ko column") as pbar:
        df_kegg_ko["KEGG_ko"] = df_kegg_ko["KEGG_ko"].str.replace(
            r"ko:(K\d+)", r"SAMPLE \1", regex=True
        )
        df_kegg_ko["KEGG_ko"] = df_kegg_ko["KEGG_ko"].str.replace(",", "\n")
        pbar.update(1)

        # Pipe text stream representations inside memory via IO StringIO buffers
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

        # Evict structural double-quote blocks artifact components
        content: str = buffer.read().replace('"', "")
        pbar.update(1)

    # Flush structural stream records into local system temp storage
    parsed_filtered_file: Path = temp_path / "parsed_KO_terms.txt"
    parsed_filtered_file.write_text(content, encoding="utf-8")

    return parsed_filtered_file


def run_kegg_decoder(
    input_file: Path | str, output_folder: Path | str, sample_name: str
) -> Path:
    """Execute KEGG-Decoder background run via atomic subprocess transactions.

    Args:
        input_file: Source path pointing to parsed KO maps text documents.
        output_folder: Export folder container directory.
        sample_name: Target text label descriptor mapping axis definitions.

    Returns:
        Path: Location pointer to the reconstructed metabolic pathway matrix table.
    """
    in_path: Path = Path(input_file)
    out_dir: Path = Path(output_folder)
    output_file: Path = out_dir / f"{sample_name}_pathways.tsv"

    # Resolve active internal package structural tree layers securely
    package_dir: Path = Path(__file__).resolve().parent
    kegg_decoder_script: Path = package_dir / "KEGG_decoder.py"
    _ensure_kegg_decoder(kegg_decoder_script)

    # Deploy structural background runtime call using subprocess environments
    with tqdm(total=1, desc="Decoding KO terms") as pbar:
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

    # Read generated output to post-process layout metadata records
    lines: list[str] = output_file.read_text(encoding="utf-8").splitlines(keepends=True)

    # Capitalize taxonomic functional metadata keys headers
    if lines:
        first_line: list[str] = lines[0].strip().split("\t")
        processed_first_line: list[str] = [
            x.capitalize() if isinstance(x, str) and x.islower() else x
            for x in first_line
        ]
        lines[0] = "\t".join(processed_first_line) + "\n"

    # Evict general placeholder indexes in favor of user sample definitions
    content: str = "".join(lines).replace("SAMPLE", sample_name)
    output_file.write_text(content, encoding="utf-8")

    return output_file
