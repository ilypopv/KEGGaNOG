"""Tests for app.py pure helpers and per-route unit tests.

Sections
--------
1. _normalize_job_id   — UUID normalisation and rejection of non-UUID strings.
2. _safe_client_filename — path-traversal defence and default-fallback logic.
3. /status route       — 404 on unknown / invalid UUID.
4. /preview route      — 404 / 400 state-machine guards.
5. /download route     — 404 / 400 state-machine guards.
6. /samples route      — 404 / 400 / 200 with real job state.
7. /pathways route     — 404 / 400 / 200 with real job state.
8. /viz route          — 422 / 404 / 400 guards + async barplot smoke.
9. Upload limits       — 413 on oversized single and multi uploads.
"""

from __future__ import annotations

import uuid
import zipfile
from pathlib import Path

import httpx2
import pytest
from fastapi.testclient import TestClient

from kegganog.app import (
    _jobs,
    _normalize_job_id,
    _safe_client_filename,
    app,
)

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def client() -> TestClient:
    return TestClient(app)


# ===========================================================================
# 1. _normalize_job_id
# ===========================================================================


class TestNormalizeJobId:
    def test_valid_uuid_is_returned_unchanged(self):
        uid = "550e8400-e29b-41d4-a716-446655440000"
        assert _normalize_job_id(uid) == uid

    def test_uuid_without_dashes_is_normalised(self):
        assert (
            _normalize_job_id("550e8400e29b41d4a716446655440000")
            == "550e8400-e29b-41d4-a716-446655440000"
        )

    def test_invalid_string_returns_none(self):
        assert _normalize_job_id("not-a-uuid") is None

    def test_empty_string_returns_none(self):
        assert _normalize_job_id("") is None


# ===========================================================================
# 2. _safe_client_filename
# ===========================================================================


class TestSafeClientFilename:
    def test_normal_filename_is_returned_unchanged(self):
        assert (
            _safe_client_filename("sample.emapper.annotations")
            == "sample.emapper.annotations"
        )

    def test_none_returns_default(self):
        assert _safe_client_filename(None) == "sample.annotations"

    def test_empty_string_returns_default(self):
        assert _safe_client_filename("") == "sample.annotations"

    def test_dot_only_returns_default(self):
        assert _safe_client_filename(".") == "sample.annotations"

    def test_dotdot_returns_default(self):
        assert _safe_client_filename("..") == "sample.annotations"

    def test_path_traversal_stripped(self):
        result = _safe_client_filename("../../etc/passwd")
        assert "/" not in result
        assert "\\" not in result
        assert result == "passwd"

    def test_leading_and_trailing_whitespace_stripped(self):
        assert _safe_client_filename("  sample.tsv  ") == "sample.tsv"


# ===========================================================================
# 3. /status route
# ===========================================================================


def test_status_unknown_uuid_returns_404(client):
    r = client.get("/status/00000000-0000-0000-0000-000000000000")
    assert r.status_code == 404


def test_status_invalid_uuid_returns_404(client):
    r = client.get("/status/not-a-uuid")
    assert r.status_code == 404


# ===========================================================================
# 4. /preview route
# ===========================================================================


def test_preview_unknown_uuid_returns_404(client):
    r = client.get("/preview/00000000-0000-0000-0000-000000000000")
    assert r.status_code == 404


def test_preview_invalid_uuid_returns_404(client):
    r = client.get("/preview/invalid-uuid")
    assert r.status_code == 404


def test_preview_job_not_done_returns_400(running_job):
    job_id = running_job
    with TestClient(app) as c:
        r = c.get(f"/preview/{job_id}")
    assert r.status_code == 400


# ===========================================================================
# 5. /download route
# ===========================================================================


def test_download_unknown_uuid_returns_404(client):
    r = client.get("/download/00000000-0000-0000-0000-000000000000")
    assert r.status_code == 404


def test_download_invalid_uuid_returns_404(client):
    r = client.get("/download/not-a-uuid")
    assert r.status_code == 404


def test_download_job_not_done_returns_400(pending_job):
    job_id = pending_job
    with TestClient(app) as c:
        r = c.get(f"/download/{job_id}")
    assert r.status_code == 400


# ===========================================================================
# 6. /samples route
# ===========================================================================


def test_samples_unknown_uuid_returns_404(client):
    r = client.get("/samples/00000000-0000-0000-0000-000000000000")
    assert r.status_code == 404


def test_samples_invalid_uuid_returns_404(client):
    r = client.get("/samples/not-a-uuid")
    assert r.status_code == 404


def test_samples_job_not_done_returns_400(running_job):
    job_id = running_job
    # Override samples to confirm it's the status guard, not missing data
    _jobs[job_id]["samples"] = ["S1"]
    with TestClient(app) as c:
        r = c.get(f"/samples/{job_id}")
    assert r.status_code == 400


def test_samples_done_returns_sample_list(registered_job):
    job_id, job = registered_job
    job["samples"] = ["Sample1", "Sample2"]
    with TestClient(app) as c:
        r = c.get(f"/samples/{job_id}")
    assert r.status_code == 200
    assert r.json()["samples"] == ["Sample1", "Sample2"]


# ===========================================================================
# 7. /pathways route
# ===========================================================================


def test_pathways_unknown_uuid_returns_404(client):
    r = client.get("/pathways/00000000-0000-0000-0000-000000000000")
    assert r.status_code == 404


def test_pathways_invalid_uuid_returns_404(client):
    r = client.get("/pathways/not-a-uuid")
    assert r.status_code == 404


def test_pathways_job_not_done_returns_400(pending_job):
    job_id = pending_job
    _jobs[job_id]["pathways"] = ["p1"]
    with TestClient(app) as c:
        r = c.get(f"/pathways/{job_id}")
    assert r.status_code == 400


def test_pathways_done_returns_pathway_list(registered_job):
    job_id, job = registered_job
    job["pathways"] = ["glycolysis", "TCA Cycle"]
    with TestClient(app) as c:
        r = c.get(f"/pathways/{job_id}")
    assert r.status_code == 200
    assert "glycolysis" in r.json()["pathways"]


# ===========================================================================
# 8. /viz route
# ===========================================================================


def test_viz_invalid_uuid_returns_404(client):
    r = client.post("/viz/not-a-uuid", data={"plot_type": "barplot"})
    assert r.status_code == 404


def test_viz_unknown_uuid_returns_404(client):
    r = client.post(
        "/viz/00000000-0000-0000-0000-000000000000", data={"plot_type": "barplot"}
    )
    assert r.status_code == 404


def test_viz_invalid_plot_type_returns_422(registered_job, client):
    job_id, _ = registered_job
    with TestClient(app) as c:
        r = c.post(f"/viz/{job_id}", data={"plot_type": "invalid_type"})
    assert r.status_code == 422


def test_viz_job_not_done_returns_400(running_job):
    job_id = running_job
    with TestClient(app) as c:
        r = c.post(f"/viz/{job_id}", data={"plot_type": "barplot"})
    assert r.status_code == 400


def test_viz_no_tsv_path_returns_400(registered_job):
    job_id, job = registered_job
    job["tsv_path"] = None
    with TestClient(app) as c:
        r = c.post(f"/viz/{job_id}", data={"plot_type": "barplot"})
    assert r.status_code == 400


@pytest.mark.asyncio
async def test_viz_barplot_with_fixture_tsv(tmp_path):
    """Barplot viz with a real TSV fixture returns a PNG response."""
    tsv_path = FIXTURES / "simple_decoder.tsv"

    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")

    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["TEST"],
        "pathways": ["a1", "a2"],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(f"/viz/{job_id}", data={"plot_type": "barplot"})
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


# ===========================================================================
# 8b. /viz route — _blocking_viz branches (corrnet, radarplot, stackedbar,
#                  streamgraph, heatmap, unknown plot_type error)
# ===========================================================================


@pytest.mark.asyncio
async def test_viz_corrnet_returns_png(tmp_path):
    tsv_path = FIXTURES / "simple_decoder.tsv"
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["TEST"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(f"/viz/{job_id}", data={"plot_type": "corrnet"})
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


@pytest.mark.asyncio
async def test_viz_corrnet_invalid_edge_cmap_falls_back_gracefully(tmp_path):
    """edge_cmap that fails regex → colormaps['coolwarm'] fallback branch (line 453)."""
    tsv_path = FIXTURES / "simple_decoder.tsv"
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["TEST"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(
                f"/viz/{job_id}",
                data={"plot_type": "corrnet", "edge_cmap": "not a valid cmap name!"},
            )
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


@pytest.mark.asyncio
async def test_viz_radarplot_no_pathways_selected_returns_500(tmp_path):
    """radarplot without pathways_selected hits the ValueError branch (lines 537-538)."""
    tsv_path = FIXTURES / "simple_decoder.tsv"
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["TEST"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(
                f"/viz/{job_id}",
                data={"plot_type": "radarplot", "pathways_selected": ""},
            )
        assert r.status_code == 500
    finally:
        _jobs.pop(job_id, None)


@pytest.mark.asyncio
async def test_viz_stackedbar_returns_png(tmp_path):
    tsv_path = FIXTURES / "multi_decoder.tsv"  # ← было simple_decoder.tsv
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["TEST"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(f"/viz/{job_id}", data={"plot_type": "stackedbar"})
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


@pytest.mark.asyncio
async def test_viz_streamgraph_returns_png(tmp_path):
    tsv_path = FIXTURES / "multi_decoder.tsv"  # ← было simple_decoder.tsv
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["TEST"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(f"/viz/{job_id}", data={"plot_type": "streamgraph"})
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


@pytest.mark.asyncio
async def test_viz_heatmap_single_simple_returns_png(tmp_path):
    """Single-sample simple heatmap — is_multi=False, heatmap_group=False (line ~870)."""
    tsv_path = FIXTURES / "simple_decoder.tsv"
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["TEST"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(
                f"/viz/{job_id}",
                data={
                    "plot_type": "heatmap",
                    "heatmap_sample_name": "TEST",
                    "heatmap_group": "false",
                },
            )
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


@pytest.mark.asyncio
async def test_viz_heatmap_single_grouped_returns_png(tmp_path):
    """Single-sample grouped heatmap — is_multi=False, heatmap_group=True."""
    tsv_path = FIXTURES / "grouped_decoder.tsv"
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["TEST"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(
                f"/viz/{job_id}",
                data={
                    "plot_type": "heatmap",
                    "heatmap_sample_name": "TEST",
                    "heatmap_group": "true",
                },
            )
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


@pytest.mark.asyncio
async def test_viz_heatmap_multi_simple_returns_png(tmp_path):
    """Multi-sample simple heatmap — is_multi=True (Function column present), group=False."""
    tsv_path = FIXTURES / "multi_decoder.tsv"
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["S1", "S2"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(
                f"/viz/{job_id}",
                data={"plot_type": "heatmap", "heatmap_group": "false"},
            )
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


@pytest.mark.asyncio
async def test_viz_heatmap_multi_grouped_returns_png(tmp_path):
    """Multi-sample grouped heatmap — is_multi=True, group=True."""
    tsv_path = FIXTURES / "multi_decoder.tsv"
    job_id = str(uuid.uuid4())
    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("dummy.txt", "x")
    _jobs[job_id] = {
        "status": "done",
        "path": str(zip_path),
        "png_path": str(tmp_path / "preview.png"),
        "tsv_path": str(tsv_path),
        "samples": ["S1", "S2"],
        "pathways": [],
        "message": "",
    }
    try:
        transport = httpx2.ASGITransport(app=app)
        async with httpx2.AsyncClient(
            transport=transport, base_url="http://test"
        ) as ac:
            r = await ac.post(
                f"/viz/{job_id}",
                data={"plot_type": "heatmap", "heatmap_group": "true"},
            )
        assert r.status_code == 200
        assert r.headers.get("content-type", "").startswith("image/png")
    finally:
        _jobs.pop(job_id, None)


# ===========================================================================
# 9. Upload limits
# ===========================================================================


def test_single_upload_oversized_returns_413(client):
    from kegganog.app import MAX_UPLOAD_BYTES_SINGLE

    oversized = b"x" * (MAX_UPLOAD_BYTES_SINGLE + 1)
    r = client.post(
        "/run",
        files={"file": ("big.tsv", oversized, "text/tab-separated-values")},
        data={"dpi": 300},
    )
    assert r.status_code == 413


def test_multi_upload_too_many_files_returns_413(client):
    from kegganog.app import MAX_MULTI_UPLOAD_FILES

    files = [
        ("files", (f"f{i}.tsv", b"x", "text/tab-separated-values"))
        for i in range(MAX_MULTI_UPLOAD_FILES + 1)
    ]
    r = client.post("/run-multi", files=files, data={})
    assert r.status_code == 413


# ===========================================================================
# 9b. /run-multi — ValidationError and batch size limit branches
# ===========================================================================


def test_run_multi_invalid_color_returns_422(client):
    """ValidationError branch in /run-multi (lines 196-205 equivalent for multi)."""
    r = client.post(
        "/run-multi",
        files=[("files", ("f.tsv", b"x", "text/tab-separated-values"))],
        data={"color": "NotAColor"},
    )
    assert r.status_code == 422


def test_run_multi_batch_exceeds_total_limit_returns_413(client):
    """Batch total > MAX_MULTI_BATCH_BYTES triggers 413 (lines 222-226)."""
    from kegganog.app import MAX_MULTI_BATCH_BYTES

    big_chunk = b"x" * (MAX_MULTI_BATCH_BYTES + 1)
    r = client.post(
        "/run-multi",
        files=[("files", ("big.tsv", big_chunk, "text/tab-separated-values"))],
        data={},
    )
    assert r.status_code == 413


def test_run_single_invalid_color_returns_422(client):
    """ValidationError branch in /run (lines 196-200)."""
    r = client.post(
        "/run",
        files={"file": ("f.tsv", b"x", "text/tab-separated-values")},
        data={"color": "NotAColor", "dpi": 300},
    )
    assert r.status_code == 422
