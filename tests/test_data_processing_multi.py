"""Tests for kegganog.processing.data_processing_multi (multi-sample path).

Sections
--------
1. parse_emapper    — prefix formatting, output filename, dash-skip, missing-column error.
2. run_kegg_decoder — sample-prefix substitution, header capitalisation, path contract.
3. merge_outputs    — empty dir, single sample, two samples, column contract.
"""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from kegganog.processing import data_processing_multi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_EMAPPER_HEADER = "## header line 1\n## header line 2\n## header line 3\n## header line 4\nKEGG_ko\tOther"


def _write_emapper(path: Path, ko_entries: list[str]) -> None:
    path.write_text(_EMAPPER_HEADER + "\n" + "\n".join(ko_entries) + "\n")


def _write_pathways_tsv(path: Path, prefix: str, values: list[float]) -> None:
    """Write a minimal KEGG-Decoder-style pathways TSV for merge tests."""
    headers = "Glycolysis\tTCA Cycle"
    vals = "\t".join(str(v) for v in values)
    path.write_text(f"{prefix}\t{headers}\n{prefix}\t{vals}\n")


# ===========================================================================
# 1. parse_emapper
# ===========================================================================


def test_parse_emapper_multi_formats_lines_with_prefix(tmp_path):
    input_file = tmp_path / "sample.tsv"
    sample_folder = tmp_path / "sample_folder"
    sample_folder.mkdir()
    _write_emapper(input_file, ["ko:K00001,ko:K00002\tX", "-\tY", "ko:K00003\tZ"])

    output_path = data_processing_multi.parse_emapper(
        str(input_file), str(sample_folder), "Sample1"
    )

    lines = Path(output_path).read_text().strip().splitlines()
    assert lines == ["Sample1 K00001", "Sample1 K00002", "Sample1 K00003"]


def test_parse_emapper_multi_output_filename_contains_prefix(tmp_path):
    input_file = tmp_path / "s.tsv"
    sample_folder = tmp_path / "folder"
    sample_folder.mkdir()
    _write_emapper(input_file, ["ko:K00001\tX"])

    output_path = data_processing_multi.parse_emapper(
        str(input_file), str(sample_folder), "SampleABC"
    )

    assert os.path.basename(output_path) == "SampleABC_parsed_KO_terms.txt"


def test_parse_emapper_multi_skips_all_dash_entries(tmp_path):
    input_file = tmp_path / "s.tsv"
    sample_folder = tmp_path / "folder"
    sample_folder.mkdir()
    _write_emapper(input_file, ["-\tX", "-\tY"])

    output_path = data_processing_multi.parse_emapper(
        str(input_file), str(sample_folder), "S1"
    )

    assert Path(output_path).read_text().strip() == ""


def test_parse_emapper_multi_missing_kegg_ko_column_raises(tmp_path):
    input_file = tmp_path / "bad.tsv"
    sample_folder = tmp_path / "folder"
    sample_folder.mkdir()
    input_file.write_text("## header\nOtherCol\tAnotherCol\nval1\tval2\n")

    with pytest.raises(KeyError, match="KEGG_ko"):
        data_processing_multi.parse_emapper(str(input_file), str(sample_folder), "S1")


# ===========================================================================
# 2. run_kegg_decoder (multi)
# ===========================================================================


def test_run_kegg_decoder_multi_substitutes_prefix(tmp_path):
    input_file = tmp_path / "ko_terms.txt"
    input_file.write_text("Sample1 K00001\n")
    file_prefix = "Sample1"
    _write_pathways_tsv(
        tmp_path / f"{file_prefix}_pathways.tsv", file_prefix, [0.5, 0.7]
    )

    with (
        patch("kegganog.processing.data_processing_multi._ensure_kegg_decoder"),
        patch("subprocess.run"),
    ):
        result = data_processing_multi.run_kegg_decoder(
            str(input_file), str(tmp_path), file_prefix
        )

    content = Path(result).read_text()
    assert "Sample1" in content
    assert "SAMPLE" not in content


def test_run_kegg_decoder_multi_capitalizes_lowercase_headers(tmp_path):
    input_file = tmp_path / "ko_terms.txt"
    input_file.write_text("S2 K00001\n")
    file_prefix = "S2"
    (tmp_path / f"{file_prefix}_pathways.tsv").write_text(
        "SAMPLE\tglycoly\tTCA Cycle\nSAMPLE\t0.4\t0.9\n"
    )

    with (
        patch("kegganog.processing.data_processing_multi._ensure_kegg_decoder"),
        patch("subprocess.run"),
    ):
        result = data_processing_multi.run_kegg_decoder(
            str(input_file), str(tmp_path), file_prefix
        )

    content = Path(result).read_text()
    assert "Glycoly" in content  # all-lowercase → capitalize()
    assert "TCA Cycle" in content  # mixed-case → unchanged


def test_run_kegg_decoder_multi_returns_correct_path(tmp_path):
    input_file = tmp_path / "ko.txt"
    input_file.write_text("X K00001\n")
    file_prefix = "MyPrefix"
    expected = tmp_path / f"{file_prefix}_pathways.tsv"
    expected.write_text("SAMPLE\tpathway\nSAMPLE\t1.0\n")

    with (
        patch("kegganog.processing.data_processing_multi._ensure_kegg_decoder"),
        patch("subprocess.run"),
    ):
        result = data_processing_multi.run_kegg_decoder(
            str(input_file), str(tmp_path), file_prefix
        )

    assert result == expected


# ===========================================================================
# 3. merge_outputs
# ===========================================================================


def test_merge_outputs_empty_directory_returns_empty_dataframe(tmp_path):
    result = data_processing_multi.merge_outputs(str(tmp_path))
    assert result.empty


def test_merge_outputs_single_sample_has_correct_columns(tmp_path):
    sample_dir = tmp_path / "temp_files" / "OnlySample"
    sample_dir.mkdir(parents=True)
    _write_pathways_tsv(
        sample_dir / "OnlySample_pathways.tsv", "OnlySample", [1.0, 0.0]
    )

    result = data_processing_multi.merge_outputs(str(tmp_path))

    assert "OnlySample" in result.columns
    assert len(result) == 2


def test_merge_outputs_two_samples_creates_merged_file(tmp_path):
    temp_files_dir = tmp_path / "temp_files"
    for prefix in ["S1", "S2"]:
        d = temp_files_dir / prefix
        d.mkdir(parents=True)
        _write_pathways_tsv(d / f"{prefix}_pathways.tsv", prefix, [0.5, 0.7])

    data_processing_multi.merge_outputs(str(tmp_path))

    assert (tmp_path / "merged_pathways.tsv").exists()


def test_merge_outputs_two_samples_column_and_function_contract(tmp_path):
    temp_files_dir = tmp_path / "temp_files"

    for prefix, vals in [("S1", [0.5, 0.7]), ("S2", [0.3, 0.8])]:
        d = temp_files_dir / prefix
        d.mkdir(parents=True)
        _write_pathways_tsv(d / f"{prefix}_pathways.tsv", prefix, vals)

    result = data_processing_multi.merge_outputs(str(tmp_path))

    assert "Function" in result.columns
    assert "S1" in result.columns
    assert "S2" in result.columns
    assert "Glycolysis" in result["Function"].values
    assert "TCA Cycle" in result["Function"].values
