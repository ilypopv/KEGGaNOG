"""Shared pytest configuration and fixtures for the KEGGaNOG test suite.

File layout
-----------
conftest.py                    ← shared fixtures, helpers, autouse setup  (all tests)
test_schemas.py                ← Pydantic model validation (pure, no I/O)
test_kegganog_cli.py           ← Typer CLI smoke + error-path + MultiSampleRunner tests
test_data_processing.py        ← single-sample data-processing pipeline steps
test_data_processing_multi.py  ← multi-sample data-processing pipeline steps
test_pipeline.py               ← pipeline.py orchestration (_pack, run_single, run_multi)
test_heatmaps_characterization.py ← cheatmap rendering reproducibility (PNG hash locks)
test_kgnplot.py                ← kgnplot public API smoke + parameter-validation tests
test_app_helpers.py            ← app.py pure helpers (_normalize_job_id, _safe_client_filename)
                                   and per-route unit tests that need _jobs manipulation
test_web_routes.py             ← HTTP-level FastAPI surface tests (TestClient / httpx2)
test_web_launch.py             ← kegganog.web launch + browser-open tests
"""

from __future__ import annotations

import io
import uuid
import zipfile
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pytest
from typer.testing import CliRunner

# Headless, deterministic rendering for tests — matches the web app Agg backend.
matplotlib.use("Agg")


# ---------------------------------------------------------------------------
# Autouse: silence side effects that would pollute every test
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _silence_heatmap_side_effects(monkeypatch):
    """Suppress plt.show() calls and tqdm progress output in all tests."""
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: None)

    import tqdm as tqdm_mod

    class _SilentTqdm(tqdm_mod.tqdm):
        def __init__(self, *args, **kwargs):
            kwargs["disable"] = True
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(tqdm_mod, "tqdm", _SilentTqdm)


# ---------------------------------------------------------------------------
# CLI runner
# ---------------------------------------------------------------------------


@pytest.fixture
def cli_runner() -> CliRunner:
    """Typer CliRunner shared by all CLI smoke tests."""
    return CliRunner()


# ---------------------------------------------------------------------------
# Minimal PNG bytes
# ---------------------------------------------------------------------------


def make_minimal_png() -> bytes:
    """Return a valid, minimal 1×1 PNG rendered by matplotlib (Agg backend).

    Used wherever a real PNG file is needed in tests without running the full
    pipeline (preview routes, pack-results, job state helpers, etc.).
    """
    buf = io.BytesIO()
    fig, ax = plt.subplots(figsize=(1, 1))
    ax.plot([0, 1], [0, 1])
    fig.savefig(buf, format="png", dpi=72)
    plt.close(fig)
    return buf.getvalue()


@pytest.fixture
def minimal_png_bytes() -> bytes:
    """Fixture wrapper around make_minimal_png() for use in test signatures."""
    return make_minimal_png()


# ---------------------------------------------------------------------------
# Completed-job dictionary helper
# ---------------------------------------------------------------------------


def make_done_job(workdir: Path) -> dict:
    """Build a minimal 'done' job-state dictionary backed by real files.

    Creates preview.png, pathways.tsv, and results.zip under *workdir* so
    that routes which serve those files (preview, download, viz) work without
    hitting the real pipeline.
    """
    png_bytes = make_minimal_png()

    png_path = workdir / "preview.png"
    png_path.write_bytes(png_bytes)

    tsv_path = workdir / "pathways.tsv"
    tsv_path.write_text("Function\tSAMPLE\nglycoly\t0.5\n", encoding="utf-8")

    zip_path = workdir / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("heatmap_figure.png", png_bytes)
        zf.writestr("pathways.tsv", tsv_path.read_text(encoding="utf-8"))

    return {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(png_path),
        "tsv_path": str(tsv_path),
        "samples": ["SAMPLE"],
        "pathways": ["glycolysis"],
        "message": "",
    }


# ---------------------------------------------------------------------------
# _jobs context manager fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def registered_job(tmp_path):
    """Insert a done job into _jobs and guarantee cleanup after the test.

    Yields (job_id, job_dict) so the test can inspect or mutate the state.
    Cleanup happens in a `finally` block, so it runs even if the test fails.

    Usage::

        def test_something(registered_job):
            job_id, job = registered_job
            r = client.get(f"/preview/{job_id}")
            assert r.status_code == 200
    """
    from kegganog.app import _jobs

    job_id = str(uuid.uuid4())
    job = make_done_job(tmp_path)
    _jobs[job_id] = job
    try:
        yield job_id, job
    finally:
        _jobs.pop(job_id, None)


@pytest.fixture
def pending_job():
    """Insert a pending (non-done) job into _jobs and guarantee cleanup."""
    from kegganog.app import _jobs

    job_id = str(uuid.uuid4())
    _jobs[job_id] = {
        "status": "pending",
        "path": None,
        "png_path": None,
        "message": "",
        "tsv_path": None,
        "samples": [],
        "pathways": [],
    }
    try:
        yield job_id
    finally:
        _jobs.pop(job_id, None)


@pytest.fixture
def running_job():
    """Insert a running (non-done) job into _jobs and guarantee cleanup."""
    from kegganog.app import _jobs

    job_id = str(uuid.uuid4())
    _jobs[job_id] = {
        "status": "running",
        "path": None,
        "png_path": None,
        "message": "",
        "tsv_path": None,
        "samples": [],
        "pathways": [],
    }
    try:
        yield job_id
    finally:
        _jobs.pop(job_id, None)
