"""HTTP-level characterization tests for the FastAPI surface.

No real pipeline runs here — the blocking worker is always patched out.
Tests are grouped by route and ordered: validation → not-found → happy path.
"""

from __future__ import annotations

import asyncio
import zipfile
from pathlib import Path
from unittest.mock import patch

import httpx2
import pytest
from conftest import make_minimal_png
from fastapi.testclient import TestClient

from kegganog.app import app
from kegganog.processing.pipeline import PipelineResult

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def client() -> TestClient:
    return TestClient(app)


def _fake_blocking_result(tmp_path: Path) -> dict:
    """Build a fake run_single return value backed by real files in tmp_path."""
    png_bytes = make_minimal_png()

    png_path = tmp_path / "preview.png"
    png_path.write_bytes(png_bytes)

    tsv_path = tmp_path / "pathways.tsv"
    tsv_path.write_text("X\tSAMPLE\na\t1.0\n", encoding="utf-8")

    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("heatmap_figure.png", png_bytes)
        zf.writestr("sample_pathways.tsv", tsv_path.read_text(encoding="utf-8"))

    return {
        "path": str(zip_path),
        "png_path": str(png_path),
        "tsv_path": str(tsv_path),
        "samples": ["SAMPLE"],
        "pathways": ["a"],
    }


# ===========================================================================
# GET /
# ===========================================================================


def test_index_returns_html(client: TestClient) -> None:
    r = client.get("/")
    assert r.status_code == 200
    assert "text/html" in r.headers.get("content-type", "")
    assert "KEGGaNOG" in r.text or "kegganog" in r.text.lower()


# ===========================================================================
# POST /run — validation and async happy path
# ===========================================================================


def test_run_validation_error_returns_422(client: TestClient) -> None:
    r = client.post(
        "/run",
        files={"file": ("test.tsv", b"x", "text/tab-separated-values")},
        data={"dpi": 10},  # dpi below minimum  # ty:ignore[invalid-argument-type]
    )
    assert r.status_code == 422


@pytest.mark.asyncio
async def test_run_completes_with_mocked_pipeline(tmp_path: Path) -> None:
    """Full async /run → poll /status → /preview → /download flow with patched worker."""

    fake = _fake_blocking_result(tmp_path)

    def fake_run_single(file_bytes, sample_name, dpi, color, group, output_dir=None):
        return PipelineResult(
            zip_path=fake["path"],
            png_path=fake["png_path"],
            tsv_path=fake["tsv_path"],
            samples=fake["samples"],
            pathways=fake["pathways"],
        )

    transport = httpx2.ASGITransport(app=app)
    # The patch must wrap the *entire* flow — POST + polling + assertions —
    # because the background worker runs in a thread pool and may call
    # run_single after the POST coroutine has already returned.  Exiting the
    # patch context before the worker finishes causes the real run_single to
    # be called with b"dummy" bytes, which blows up on missing KEGG_ko column.
    with patch("kegganog.app.run_single", fake_run_single):
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(
                "/run",
                files={"file": ("in.tsv", b"dummy", "text/tab-separated-values")},
                data={
                    "dpi": 300,
                    "color": "Blues",
                    "sample_name": "SAMPLE",
                    "group": "false",
                },
            )
            assert r.status_code == 200, r.text
            job_id = r.json()["job_id"]

            # Poll until done (max 5 seconds)
            status_payload = None
            for _ in range(100):
                await asyncio.sleep(0.05)
                s = await ac.get(f"/status/{job_id}")
                status_payload = s.json()
                if status_payload.get("status") in ("done", "error"):
                    break

            assert status_payload is not None
            assert status_payload["status"] == "done", status_payload

            # Preview
            prev = await ac.get(f"/preview/{job_id}")
            assert prev.status_code == 200
            assert prev.headers.get("content-type", "").startswith("image/png")

            # Download
            dl = await ac.get(f"/download/{job_id}")
            assert dl.status_code == 200
            assert "zip" in dl.headers.get("content-type", "").lower()


# ===========================================================================
# GET /status
# ===========================================================================


def test_status_not_found_returns_404(client: TestClient) -> None:
    r = client.get("/status/00000000-0000-0000-0000-000000000000")
    assert r.status_code == 404


# ===========================================================================
# POST /run-multi
# ===========================================================================


def test_run_multi_rejects_empty_files_returns_422(client: TestClient) -> None:
    r = client.post("/run-multi", files=[], data={})
    assert r.status_code == 422
