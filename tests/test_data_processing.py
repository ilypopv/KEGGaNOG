"""Tests for kegganog.processing.data_processing (single-sample path).

Sections
--------
1. parse_emapper         — KO filtering, formatting, output path.
2. _ensure_kegg_decoder  — download-on-miss, skip-if-present, error handling.
3. run_kegg_decoder      — sample-name substitution, header capitalisation, path contract.
"""

from __future__ import annotations

import pathlib
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from kegganog.processing import data_processing
from kegganog.processing.data_processing import _ensure_kegg_decoder

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_EMAPPER_HEADER = "## header line 1\n## header line 2\n## header line 3\n## header line 4\nKEGG_ko\tOther"


def _write_emapper(path: Path, ko_entries: list[str]) -> None:
    path.write_text(_EMAPPER_HEADER + "\n" + "\n".join(ko_entries) + "\n")


# ===========================================================================
# 1. parse_emapper
# ===========================================================================


def test_parse_emapper_filters_dashes_and_formats_ko_terms(tmp_path):
    input_file = tmp_path / "emapper.tsv"
    temp_folder = tmp_path / "temp"
    temp_folder.mkdir()
    _write_emapper(input_file, ["ko:K00001,ko:K00002\tX", "-\tY", "ko:K00003\tZ"])

    output_path = data_processing.parse_emapper(str(input_file), str(temp_folder))

    lines = pathlib.Path(output_path).read_text().strip().splitlines()
    assert lines == ["SAMPLE K00001", "SAMPLE K00002", "SAMPLE K00003"]


def test_parse_emapper_all_dash_entries_produces_empty_output(tmp_path):
    input_file = tmp_path / "emapper.tsv"
    temp_folder = tmp_path / "temp"
    temp_folder.mkdir()
    _write_emapper(input_file, ["-\tX", "-\tY"])

    output_path = data_processing.parse_emapper(str(input_file), str(temp_folder))

    assert pathlib.Path(output_path).read_text().strip() == ""


# ===========================================================================
# 2. _ensure_kegg_decoder
# ===========================================================================


def test_ensure_kegg_decoder_skips_download_when_file_exists(tmp_path):
    script_path = tmp_path / "KEGG_decoder.py"
    script_path.write_text("# fake")

    with patch("requests.get") as mock_get:
        _ensure_kegg_decoder(script_path)

    mock_get.assert_not_called()
    assert script_path.exists()


def test_ensure_kegg_decoder_downloads_when_missing(tmp_path):
    script_path = tmp_path / "KEGG_decoder.py"
    mock_response = MagicMock()
    mock_response.content = b"# fake kegg decoder"
    mock_response.raise_for_status = MagicMock()

    with patch("requests.get", return_value=mock_response) as mock_get:
        _ensure_kegg_decoder(script_path)

    mock_get.assert_called_once()
    assert script_path.read_bytes() == b"# fake kegg decoder"


def test_ensure_kegg_decoder_raises_on_network_failure(tmp_path):
    script_path = tmp_path / "KEGG_decoder.py"

    with (
        patch("requests.get", side_effect=Exception("network error")),
        pytest.raises(RuntimeError, match="Failed to download"),
    ):
        _ensure_kegg_decoder(script_path)

    assert not script_path.exists()


# ===========================================================================
# 3. run_kegg_decoder
# ===========================================================================


def _write_pathways_tsv(
    path: Path, sample_name: str, headers: list[str], values: list[float]
) -> None:
    path.write_text(
        "SAMPLE\t"
        + "\t".join(headers)
        + "\n"
        + "SAMPLE\t"
        + "\t".join(str(v) for v in values)
        + "\n"
    )


def test_run_kegg_decoder_substitutes_sample_name(tmp_path):
    input_file = tmp_path / "ko_terms.txt"
    input_file.write_text("SAMPLE K00001\n")
    sample_name = "MySample"
    _write_pathways_tsv(
        tmp_path / f"{sample_name}_pathways.tsv",
        sample_name,
        ["glycoly", "TCA Cycle"],
        [0.5, 0.7],
    )

    with (
        patch("kegganog.processing.data_processing._ensure_kegg_decoder"),
        patch("subprocess.run"),
    ):
        result = data_processing.run_kegg_decoder(
            str(input_file), str(tmp_path), sample_name
        )

    content = Path(result).read_text()
    assert "MySample" in content
    assert "SAMPLE" not in content


def test_run_kegg_decoder_capitalizes_all_lowercase_headers(tmp_path):
    input_file = tmp_path / "ko_terms.txt"
    input_file.write_text("SAMPLE K00001\n")
    sample_name = "S1"
    _write_pathways_tsv(
        tmp_path / f"{sample_name}_pathways.tsv",
        sample_name,
        ["glycoly", "TCA Cycle"],
        [0.5, 0.7],
    )

    with (
        patch("kegganog.processing.data_processing._ensure_kegg_decoder"),
        patch("subprocess.run"),
    ):
        result = data_processing.run_kegg_decoder(
            str(input_file), str(tmp_path), sample_name
        )

    content = Path(result).read_text()
    assert "Glycoly" in content  # all-lowercase → capitalize()
    assert "TCA Cycle" in content  # mixed-case → unchanged


def test_run_kegg_decoder_returns_correct_output_path(tmp_path):
    input_file = tmp_path / "ko_terms.txt"
    input_file.write_text("SAMPLE K00001\n")
    sample_name = "TestSample"
    expected = tmp_path / f"{sample_name}_pathways.tsv"
    _write_pathways_tsv(expected, sample_name, ["pathway"], [1.0])

    with (
        patch("kegganog.processing.data_processing._ensure_kegg_decoder"),
        patch("subprocess.run"),
    ):
        result = data_processing.run_kegg_decoder(
            str(input_file), str(tmp_path), sample_name
        )

    assert result == expected
