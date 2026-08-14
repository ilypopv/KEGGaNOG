"""Tests for kegganog.processing.pipeline.

Sections
--------
1. _secure_temp_path  — existence, suffix, uniqueness.
2. _pack_results      — zip creation, PNG presence, contents, missing-PNG error.
3. PipelineResult     — field contract and default lists.
4. run_single         — web mode (output_dir=None) and CLI mode (output_dir provided).
5. run_multi          — web mode and CLI mode.
6. _run_single_in_dir — simple and grouped heatmap branches.
7. _run_multi_in_dir  — simple and grouped heatmap branches.
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import pytest

from kegganog.processing.pipeline import (
    PipelineResult,
    _pack_results,
    _secure_temp_path,
    run_multi,
    run_single,
)

# ---------------------------------------------------------------------------
# Module-level helpers (no inline imports, no repeated matplotlib.use calls)
# ---------------------------------------------------------------------------


def _make_minimal_png_bytes() -> bytes:
    buf = io.BytesIO()
    fig, _ax = plt.subplots(figsize=(1, 1))
    fig.savefig(buf, format="png", dpi=72)
    plt.close(fig)
    return buf.getvalue()


def _make_output_dir_with_png(tmp_path: Path) -> Path:
    """Return an output dir that satisfies _pack_results (has heatmap_figure.png)."""
    out = tmp_path / "output"
    out.mkdir()
    (out / "heatmap_figure.png").write_bytes(_make_minimal_png_bytes())
    return out


def _fake_heatmap(tsv_or_df, output_folder, dpi, color, sample_name=None, **kwargs):
    """Drop-in replacement for any generate_heatmap* call; writes a real PNG."""
    (Path(output_folder) / "heatmap_figure.png").write_bytes(_make_minimal_png_bytes())


def _fake_run_multi_in_dir(named_files, dpi, color, group, output_dir):
    merged_tsv = output_dir / "merged_pathways.tsv"
    merged_tsv.write_text("Function\tS1\tS2\nglycoly\t0.5\t0.3\n")
    (output_dir / "heatmap_figure.png").write_bytes(_make_minimal_png_bytes())
    return ["S1", "S2"], ["glycoly"], str(merged_tsv)


# ===========================================================================
# 1. _secure_temp_path
# ===========================================================================


def test_secure_temp_path_has_correct_suffix():
    p = _secure_temp_path(".tsv")
    assert p.suffix == ".tsv"
    p.unlink()


def test_secure_temp_path_creates_file():
    p = _secure_temp_path(".tsv")
    assert p.exists()
    p.unlink()


def test_secure_temp_path_produces_unique_paths():
    p1, p2 = _secure_temp_path(".png"), _secure_temp_path(".png")
    assert p1 != p2
    p1.unlink()
    p2.unlink()


# ===========================================================================
# 2. _pack_results
# ===========================================================================


def test_pack_results_creates_zip_and_png(tmp_path):
    out = _make_output_dir_with_png(tmp_path)
    zip_path, png_path = _pack_results(out)

    assert Path(zip_path).exists()
    assert Path(png_path).exists()


def test_pack_results_output_has_correct_suffixes(tmp_path):
    out = _make_output_dir_with_png(tmp_path)
    zip_path, png_path = _pack_results(out)

    assert Path(zip_path).suffix == ".zip"
    assert Path(png_path).suffix == ".png"


def test_pack_results_zip_contains_heatmap_and_extra_files(tmp_path):
    out = _make_output_dir_with_png(tmp_path)
    (out / "extra.tsv").write_text("data")
    zip_path, _ = _pack_results(out)

    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
    assert any("heatmap_figure.png" in n for n in names)
    assert any("extra.tsv" in n for n in names)


def test_pack_results_raises_when_png_missing(tmp_path):
    out = tmp_path / "empty"
    out.mkdir()
    with pytest.raises(FileNotFoundError):
        _pack_results(out)


# ===========================================================================
# 3. PipelineResult
# ===========================================================================


def test_pipeline_result_stores_all_fields():
    r = PipelineResult(
        zip_path=Path("/tmp/a.zip"),
        png_path=Path("/tmp/a.png"),
        tsv_path=Path("/tmp/a.tsv"),
        samples=["S1", "S2"],
        pathways=["glycolysis"],
    )
    assert r.zip_path == Path("/tmp/a.zip")
    assert r.samples == ["S1", "S2"]
    assert r.pathways == ["glycolysis"]


def test_pipeline_result_default_lists_are_empty():
    r = PipelineResult(zip_path=Path("z"), png_path=Path("p"), tsv_path=Path("t"))
    assert r.samples == []
    assert r.pathways == []


# ===========================================================================
# 4. run_single
# ===========================================================================


def _make_fake_run_single_in_dir(output_dir_ref: dict):
    """Return a patched _run_single_in_dir that writes required files and records output_dir."""

    def fake(file_bytes, sample_name, dpi, color, group, output_dir):
        output_dir_ref["captured"] = output_dir
        tsv = output_dir / "SAMPLE_pathways.tsv"
        tsv.write_text("Function\ta1\na1\t0.5\n")
        (output_dir / "heatmap_figure.png").write_bytes(_make_minimal_png_bytes())
        return ["SAMPLE"], ["a1"], str(tsv)

    return fake


def test_run_single_web_mode_returns_pipeline_result():
    ref = {}
    with patch(
        "kegganog.processing.pipeline._run_single_in_dir",
        _make_fake_run_single_in_dir(ref),
    ):
        result = run_single(
            file_bytes=b"dummy",
            sample_name="SAMPLE",
            dpi=72,
            color="Blues",
            group=False,
        )

    assert isinstance(result, PipelineResult)
    assert result.samples == ["SAMPLE"]
    assert result.pathways == ["a1"]
    assert Path(result.zip_path).exists()
    assert Path(result.png_path).exists()
    assert Path(result.tsv_path).exists()


def test_run_single_cli_mode_writes_into_provided_output_dir(tmp_path):
    out = tmp_path / "output"
    out.mkdir()
    ref = {}

    with patch(
        "kegganog.processing.pipeline._run_single_in_dir",
        _make_fake_run_single_in_dir(ref),
    ):
        result = run_single(
            file_bytes=b"dummy",
            sample_name="SAMPLE",
            dpi=72,
            color="Blues",
            group=False,
            output_dir=str(out),
        )

    assert isinstance(result, PipelineResult)
    assert (out / "heatmap_figure.png").exists()


# ===========================================================================
# 5. run_multi
# ===========================================================================


def test_run_multi_web_mode_returns_pipeline_result():
    with patch(
        "kegganog.processing.pipeline._run_multi_in_dir", _fake_run_multi_in_dir
    ):
        result = run_multi(
            named_files=[("s1.annotations", b"dummy"), ("s2.annotations", b"dummy")],
            dpi=72,
            color="Blues",
            group=False,
        )

    assert isinstance(result, PipelineResult)
    assert result.samples == ["S1", "S2"]
    assert result.pathways == ["glycoly"]
    assert Path(result.zip_path).exists()


def test_run_multi_cli_mode_writes_into_provided_output_dir(tmp_path):
    out = tmp_path / "output"
    out.mkdir()

    with patch(
        "kegganog.processing.pipeline._run_multi_in_dir", _fake_run_multi_in_dir
    ):
        result = run_multi(
            named_files=[("s1.annotations", b"dummy")],
            dpi=72,
            color="Blues",
            group=False,
            output_dir=str(out),
        )

    assert isinstance(result, PipelineResult)
    assert (out / "heatmap_figure.png").exists()


# ===========================================================================
# 6. _run_single_in_dir
# ===========================================================================


def test_run_single_in_dir_simple_heatmap(tmp_path):
    from kegganog.processing.pipeline import _run_single_in_dir

    out = tmp_path / "output"
    out.mkdir()
    (out / "temp_files").mkdir()

    with (
        patch(
            "kegganog.processing.pipeline.parse_emapper",
            return_value=str(tmp_path / "p.txt"),
        ),
        patch(
            "kegganog.processing.pipeline.run_kegg_decoder",
            side_effect=lambda p, d, n: _write_fake_tsv(d, n),
        ),
        patch("kegganog.cheatmaps.simple_heatmap.generate_heatmap", _fake_heatmap),
    ):
        samples, pathways, tsv_path = _run_single_in_dir(
            file_bytes=b"dummy",
            sample_name="SAMPLE",
            dpi=72,
            color="Blues",
            group=False,
            output_dir=out,
        )

    assert samples == ["SAMPLE"]
    assert "a1" in pathways
    assert Path(tsv_path).exists()


def test_run_single_in_dir_grouped_heatmap(tmp_path):
    from kegganog.processing.pipeline import _run_single_in_dir

    out = tmp_path / "output"
    out.mkdir()

    with (
        patch(
            "kegganog.processing.pipeline.parse_emapper",
            return_value=str(tmp_path / "p.txt"),
        ),
        patch(
            "kegganog.processing.pipeline.run_kegg_decoder",
            side_effect=lambda p, d, n: _write_fake_tsv(d, n),
        ),
        patch(
            "kegganog.cheatmaps.grouped_heatmap.generate_grouped_heatmap", _fake_heatmap
        ),
    ):
        samples, _, _ = _run_single_in_dir(
            file_bytes=b"dummy",
            sample_name="SAMPLE",
            dpi=72,
            color="Blues",
            group=True,
            output_dir=out,
        )

    assert samples == ["SAMPLE"]


def _write_fake_tsv(output_folder, sample_name):
    tsv = Path(output_folder) / f"{sample_name}_pathways.tsv"
    tsv.write_text("Function\ta1\na1\t0.5\n")
    return str(tsv)


# ===========================================================================
# 7. _run_multi_in_dir
# ===========================================================================


def _make_fake_merge(output_folder):
    import pandas as pd

    df = pd.DataFrame(
        {"Function": ["glycolysis", "TCA Cycle"], "S1": [0.5, 0.7], "S2": [0.3, 0.8]}
    )
    merged = Path(output_folder) / "merged_pathways.tsv"
    df.to_csv(merged, sep="\t", index=False)
    return df


def test_run_multi_in_dir_simple_heatmap(tmp_path):
    from kegganog.processing.pipeline import _run_multi_in_dir

    out = tmp_path / "output"
    out.mkdir()

    def fake_heatmap_multi(df, output_folder, dpi, color, **kwargs):
        _fake_heatmap(df, output_folder, dpi, color)

    with (
        patch(
            "kegganog.processing.pipeline.parse_emapper_multi",
            return_value=str(tmp_path / "p.txt"),
        ),
        patch(
            "kegganog.processing.pipeline.run_kegg_decoder_multi",
            return_value=str(tmp_path / "d.tsv"),
        ),
        patch(
            "kegganog.processing.pipeline.merge_outputs", side_effect=_make_fake_merge
        ),
        patch(
            "kegganog.cheatmaps.simple_heatmap_multi.generate_heatmap_multi",
            fake_heatmap_multi,
        ),
    ):
        samples, pathways, _ = _run_multi_in_dir(
            named_files=[
                ("s1.emapper.annotations", b"x"),
                ("s2.emapper.annotations", b"x"),
            ],
            dpi=72,
            color="Blues",
            group=False,
            output_dir=out,
        )

    assert "S1" in samples
    assert "S2" in samples
    assert "glycolysis" in pathways


def test_run_multi_in_dir_grouped_heatmap(tmp_path):
    from kegganog.processing.pipeline import _run_multi_in_dir

    out = tmp_path / "output"
    out.mkdir()

    def fake_merge_single(output_folder):
        import pandas as pd

        df = pd.DataFrame({"Function": ["glycolysis"], "S1": [0.5]})
        (Path(output_folder) / "merged_pathways.tsv").write_text(
            "Function\tS1\nglycolysis\t0.5\n"
        )
        return df

    def fake_heatmap_grouped(df, output_folder, dpi, color, **kwargs):
        _fake_heatmap(df, output_folder, dpi, color)

    with (
        patch(
            "kegganog.processing.pipeline.parse_emapper_multi",
            return_value=str(tmp_path / "p.txt"),
        ),
        patch(
            "kegganog.processing.pipeline.run_kegg_decoder_multi",
            return_value=str(tmp_path / "d.tsv"),
        ),
        patch(
            "kegganog.processing.pipeline.merge_outputs", side_effect=fake_merge_single
        ),
        patch(
            "kegganog.cheatmaps.grouped_heatmap_multi.generate_grouped_heatmap_multi",
            fake_heatmap_grouped,
        ),
    ):
        samples, _, _ = _run_multi_in_dir(
            named_files=[("s1.emapper.annotations", b"x")],
            dpi=72,
            color="Blues",
            group=True,
            output_dir=out,
        )

    assert samples == ["S1"]
