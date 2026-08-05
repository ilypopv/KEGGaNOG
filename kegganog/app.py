"""FastAPI web application orchestration gateway for KEGGaNOG pipelines.

This module exposes asynchronous endpoints for single- and multi-sample functional
profile analysis, manages long-running background tasks via thread pools, and serves
customizable Matplotlib/Seaborn visualization plots.
"""

import matplotlib

# Force non-interactive Agg backend before any matplotlib sub-modules load
matplotlib.use("Agg")

import asyncio
import logging
import os
import shutil
import tempfile
import uuid
from importlib.metadata import version as _metadata_version
from pathlib import Path
from typing import Literal

import pandas as pd
from fastapi import FastAPI, Form, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from matplotlib import colormaps
from pydantic import ValidationError
from starlette.responses import Response

from .processing.pipeline import run_multi, run_single
from .schemas import JobStatus, ValidColor, WebParams

FontWeight = Literal["normal", "bold", "heavy", "light"]
FontStyle = Literal["normal", "italic", "oblique"]
SortOrder = Literal["ascending", "descending"]
LineStyle = Literal["solid", "dashed", "dashdot", "dotted", "-"]
HorizontalAlignment = Literal["left", "right", "center"]

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# FastAPI application instance
# ---------------------------------------------------------------------------

app = FastAPI(
    title="KEGGaNOG",
    version=_metadata_version("kegganog"),
    docs_url=None,
    redoc_url=None,
)

_static_dir = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=_static_dir), name="static")

# ---------------------------------------------------------------------------
# In-memory job store
# ---------------------------------------------------------------------------

# Explicit typing for central tracking registry mapping job UUIDs to properties
_jobs: dict[str, dict[str, str | Path | None | list[str]]] = {}

# Upload constraints and boundary allocations
MAX_UPLOAD_BYTES_SINGLE: int = 80 * 1024 * 1024
MAX_UPLOAD_BYTES_PER_MULTI_FILE: int = 80 * 1024 * 1024
MAX_MULTI_UPLOAD_FILES: int = 64
MAX_MULTI_BATCH_BYTES: int = 400 * 1024 * 1024

_ALLOWED_PLOT_TYPES: frozenset[str] = frozenset(
    {"barplot", "corrnet", "radarplot", "stackedbar", "streamgraph", "heatmap"}
)


def _secure_temp_path(suffix: str) -> Path:
    """Return a closed secure file path inside the OS temporary storage container."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    return Path(path)


def _safe_client_filename(filename: str | None) -> str:
    """Sanitize path-traversal anomalies out of multi-part uploaded file streams."""
    raw = (filename or "").strip() or "sample.annotations"
    base = Path(raw).name
    if not base or base in {".", ".."}:
        base = "sample.annotations"
    return base


async def _read_upload_with_limit(upload: UploadFile, max_bytes: int) -> bytes:
    """Read binary chunks sequentially while throwing exceptions if bounds are exceeded."""
    chunks: list[bytes] = []
    total = 0
    chunk_size = 1024 * 1024
    while True:
        chunk = await upload.read(chunk_size)
        if not chunk:
            break
        total += len(chunk)
        if total > max_bytes:
            raise ValueError(
                f"Uploaded file exceeds maximum size ({max_bytes // (1024 * 1024)} MiB)."
            )
        chunks.append(chunk)
    return b"".join(chunks)


def _normalize_job_id(job_id: str) -> str | None:
    """Validate input character strings against canonical UUID standards."""
    try:
        return str(uuid.UUID(job_id))
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Route 1: Serve the HTML user interface
# ---------------------------------------------------------------------------


@app.get("/", response_class=HTMLResponse)
async def index() -> HTMLResponse:
    """Render index interface dashboards injection pipeline version variables."""
    html_file = Path(__file__).parent / "static" / "index.html"
    html = html_file.read_text(encoding="utf-8").replace(
        "__VERSION__", _metadata_version("kegganog")
    )
    return HTMLResponse(html)


# ---------------------------------------------------------------------------
# Route 2: Single-sample functional profiling analysis
# ---------------------------------------------------------------------------


@app.post("/run", response_model=JobStatus)
async def run_analysis(
    file: UploadFile,
    dpi: int = Form(default=300),
    color: ValidColor = Form(default="Blues"),
    sample_name: str = Form(default="SAMPLE"),
    group: bool = Form(default=False),
) -> JobStatus | JSONResponse:
    """Validate incoming single file parameters and schedule long-running threads."""
    try:
        params = WebParams(dpi=dpi, color=color, sample_name=sample_name, group=group)
    except ValidationError as e:
        messages = [f"{err['loc'][0]}: {err['msg']}" for err in e.errors()]
        return JSONResponse(status_code=422, content={"detail": messages})

    try:
        file_bytes = await _read_upload_with_limit(file, MAX_UPLOAD_BYTES_SINGLE)
    except ValueError as e:
        return JSONResponse(status_code=413, content={"detail": [str(e)]})

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
    asyncio.create_task(_run_job(job_id, file_bytes, params))
    return JobStatus(job_id=job_id, status="pending")


# ---------------------------------------------------------------------------
# Route 3: Multi-sample batch analysis
# ---------------------------------------------------------------------------


@app.post("/run-multi", response_model=JobStatus)
async def run_analysis_multi(
    files: list[UploadFile],
    dpi: int = Form(default=300),
    color: ValidColor = Form(default="Blues"),
    group: bool = Form(default=False),
) -> JobStatus | JSONResponse:
    """Process an uploaded matrix collection array under cumulative memory safeguards."""
    valid_files = [f for f in files if f.filename and f.filename.strip()]
    if not valid_files:
        return JSONResponse(
            status_code=422,
            content={"detail": ["files: at least one valid file required."]},
        )
    if len(valid_files) > MAX_MULTI_UPLOAD_FILES:
        return JSONResponse(
            status_code=413,
            content={"detail": [f"Too many files (max {MAX_MULTI_UPLOAD_FILES})."]},
        )
    try:
        params = WebParams(dpi=dpi, color=color, sample_name="MULTI", group=group)
    except ValidationError as e:
        messages = [f"{err['loc'][0]}: {err['msg']}" for err in e.errors()]
        return JSONResponse(status_code=422, content={"detail": messages})

    named_files: list[tuple[str, bytes]] = []
    batch_total = 0
    for f in valid_files:
        try:
            data = await _read_upload_with_limit(f, MAX_UPLOAD_BYTES_PER_MULTI_FILE)
        except ValueError as e:
            return JSONResponse(status_code=413, content={"detail": [str(e)]})
        batch_total += len(data)
        if batch_total > MAX_MULTI_BATCH_BYTES:
            return JSONResponse(
                status_code=413,
                content={"detail": ["Total upload size exceeds server limit."]},
            )
        named_files.append((_safe_client_filename(f.filename), data))

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
    asyncio.create_task(_run_job_multi(job_id, named_files, params))
    return JobStatus(job_id=job_id, status="pending")


# ---------------------------------------------------------------------------
# Route 4: Poll job compilation status
# ---------------------------------------------------------------------------


@app.get("/status/{job_id}")
async def get_status(job_id: str) -> JSONResponse:
    """Return runtime step summaries execution states or failure details."""
    nid = _normalize_job_id(job_id)
    if nid is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    job = _jobs.get(nid)
    if job is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    return JSONResponse(
        content={"job_id": nid, "status": job["status"], "message": job["message"]}
    )


# ---------------------------------------------------------------------------
# Route 5: Serve the PNG canvas for inline client-side web rendering
# ---------------------------------------------------------------------------


@app.get("/preview/{job_id}")
async def preview(job_id: str) -> Response:
    """Transmit primary image streams to frontend frame nodes once flagged done."""
    nid = _normalize_job_id(job_id)
    if nid is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    job = _jobs.get(nid)
    if job is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    if job["status"] != "done":
        return JSONResponse(
            status_code=400, content={"detail": "Analysis not finished yet."}
        )
    return FileResponse(path=str(job["png_path"]), media_type="image/png")


# ---------------------------------------------------------------------------
# Route 6: Download the compressed analytical output bundle zip archive
# ---------------------------------------------------------------------------


@app.get("/download/{job_id}")
async def download(job_id: str) -> Response:
    """Pack summary archives down to pipeline consumers on demand."""
    nid = _normalize_job_id(job_id)
    if nid is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    job = _jobs.get(nid)
    if job is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    if job["status"] != "done":
        return JSONResponse(
            status_code=400, content={"detail": "Analysis not finished yet."}
        )
    return FileResponse(
        path=str(job["path"]),
        media_type="application/zip",
        filename=f"kegganog_{nid[:8]}.zip",
    )


# ---------------------------------------------------------------------------
# Route 7: Return sample structural identifiers index maps
# ---------------------------------------------------------------------------


@app.get("/samples/{job_id}")
async def get_samples(job_id: str) -> JSONResponse:
    """Expose tracked list configurations for dynamic multiselect dropdown lists."""
    nid = _normalize_job_id(job_id)
    if nid is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    job = _jobs.get(nid)
    if job is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    if job["status"] != "done":
        return JSONResponse(
            status_code=400, content={"detail": "Analysis not finished yet."}
        )
    return JSONResponse(content={"samples": job["samples"]})


# ---------------------------------------------------------------------------
# Route 8: Return parsed category pathways catalogs lists
# ---------------------------------------------------------------------------


@app.get("/pathways/{job_id}")
async def get_pathways(job_id: str) -> JSONResponse:
    """Expose extracted keys facilitating radar axis filtering matrices updates."""
    nid = _normalize_job_id(job_id)
    if nid is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    job = _jobs.get(nid)
    if job is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    if job["status"] != "done":
        return JSONResponse(
            status_code=400, content={"detail": "Analysis not finished yet."}
        )
    return JSONResponse(content={"pathways": job["pathways"]})


# ---------------------------------------------------------------------------
# Route 9: Execute reactive graphic rendering subroutines over static caches
# ---------------------------------------------------------------------------


@app.post("/viz/{job_id}")
async def run_viz(
    job_id: str,
    plot_type: str = Form(...),
    dpi: int = Form(default=300),
    figwidth: int = Form(default=0),
    figheight: int = Form(default=0),
    title: str = Form(default=""),
    title_fontsize: float = Form(default=16.0),
    title_color: str = Form(default="black"),
    title_weight: FontWeight = Form(default="normal"),
    title_style: FontStyle = Form(default="normal"),
    background_color: str = Form(default="white"),
    xlabel: str = Form(default=""),
    xlabel_fontsize: float = Form(default=14.0),
    xlabel_color: str = Form(default="black"),
    xlabel_weight: FontWeight = Form(default="normal"),
    xlabel_style: FontStyle = Form(default="normal"),
    ylabel: str = Form(default=""),
    ylabel_fontsize: float = Form(default=14.0),
    ylabel_color: str = Form(default="black"),
    ylabel_weight: FontWeight = Form(default="normal"),
    ylabel_style: FontStyle = Form(default="normal"),
    xticks_fontsize: float = Form(default=12.0),
    xticks_color: str = Form(default="black"),
    xticks_weight: FontWeight = Form(default="normal"),
    xticks_style: FontStyle = Form(default="normal"),
    xticks_rotation: float = Form(default=0.0),
    xticks_ha: HorizontalAlignment = Form(default="center"),
    yticks_fontsize: float = Form(default=12.0),
    yticks_color: str = Form(default="black"),
    yticks_weight: FontWeight = Form(default="normal"),
    yticks_style: FontStyle = Form(default="normal"),
    grid: bool = Form(default=True),
    grid_linestyle: str = Form(default="--"),
    grid_alpha: float = Form(default=0.7),
    cmap: str = Form(default=""),
    cmap_range_min: int = Form(default=8),
    cmap_range_max: int = Form(default=30),
    sort_order: SortOrder = Form(default="descending"),
    box_color: str = Form(default="blue"),
    showfliers: bool = Form(default=True),
    grid_color: str = Form(default="gray"),
    grid_linewidth: float = Form(default=0.5),
    threshold: float = Form(default=0.5),
    node_size: float = Form(default=700.0),
    node_color: str = Form(default="#A3D5FF"),
    node_edgecolors: str = Form(default="#03045E"),
    node_linewidths: float = Form(default=1.5),
    label_fontsize: float = Form(default=8.0),
    label_color: str = Form(default="#03045E"),
    label_weight: FontWeight = Form(default="normal"),
    edge_cmap: str = Form(default="coolwarm"),
    cbar_size: float = Form(default=0.5),
    pathways_selected: str = Form(default=""),
    colors_selected: str = Form(default=""),
    sample_order: str = Form(default=""),
    fill_alpha: float = Form(default=0.25),
    line_width: float = Form(default=2.0),
    line_style: LineStyle = Form(default="solid"),
    label_background: str = Form(default=""),
    label_edgecolor: str = Form(default=""),
    label_pad: float = Form(default=1.05),
    show_legend: bool = Form(default=True),
    legend_loc: str = Form(default="upper right"),
    bar_width: float = Form(default=0.6),
    edgecolor: str = Form(default="black"),
    edge_linewidth: float = Form(default=0.3),
    stream_fill_alpha: float = Form(default=1.0),
    legend_fontsize: float = Form(default=9.0),
    legend_bbox_x: float = Form(default=1.05),
    legend_bbox_y: float = Form(default=1.0),
    heatmap_color: str = Form(default="Blues"),
    heatmap_group: bool = Form(default=False),
    heatmap_sample_name: str = Form(default=""),
    heatmap_dpi: int = Form(default=300),
) -> Response:
    """Parse runtime plotting parameters and delegate drawing work to worker threads."""
    import json

    if plot_type not in _ALLOWED_PLOT_TYPES:
        return JSONResponse(
            status_code=422,
            content={"detail": [f"Invalid plot_type: {plot_type!r}."]},
        )

    nid = _normalize_job_id(job_id)
    if nid is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    job = _jobs.get(nid)
    if job is None:
        return JSONResponse(status_code=404, content={"detail": "Job not found."})
    if job["status"] != "done":
        return JSONResponse(
            status_code=400, content={"detail": "Analysis not finished yet."}
        )

    tsv_path_val = job.get("tsv_path")
    if not tsv_path_val:
        return JSONResponse(
            status_code=400, content={"detail": "No TSV data available for this job."}
        )

    def _parse_json_list(s: str) -> list:
        s_clean = s.strip()
        if not s_clean:
            return []
        try:
            return json.loads(s_clean)
        except Exception:
            return []

    figsize = (
        (float(figwidth), float(figheight)) if figwidth > 0 and figheight > 0 else None
    )
    loop = asyncio.get_event_loop()

    try:
        png_path = await loop.run_in_executor(
            None,
            _blocking_viz,
            str(tsv_path_val),
            plot_type,
            figsize,
            dpi,
            heatmap_color,
            heatmap_group,
            heatmap_sample_name,
            heatmap_dpi,
            title,
            title_fontsize,
            title_color,
            title_weight,
            title_style,
            background_color,
            xlabel,
            xlabel_fontsize,
            xlabel_color,
            xlabel_weight,
            xlabel_style,
            ylabel,
            ylabel_fontsize,
            ylabel_color,
            ylabel_weight,
            ylabel_style,
            xticks_fontsize,
            xticks_color,
            xticks_weight,
            xticks_style,
            xticks_rotation,
            xticks_ha,
            yticks_fontsize,
            yticks_color,
            yticks_weight,
            yticks_style,
            grid,
            grid_linestyle,
            grid_alpha,
            cmap,
            cmap_range_min,
            cmap_range_max,
            sort_order,
            box_color,
            showfliers,
            grid_color,
            grid_linewidth,
            threshold,
            node_size,
            node_color,
            node_edgecolors,
            node_linewidths,
            label_fontsize,
            label_color,
            label_weight,
            edge_cmap,
            cbar_size,
            _parse_json_list(pathways_selected),
            _parse_json_list(colors_selected),
            _parse_json_list(sample_order),
            fill_alpha,
            line_width,
            line_style,
            label_background or None,
            label_edgecolor or None,
            label_pad,
            show_legend,
            legend_loc,
            bar_width,
            edgecolor,
            edge_linewidth,
            stream_fill_alpha,
            legend_fontsize,
            (legend_bbox_x, legend_bbox_y),
        )
    except Exception as e:
        return JSONResponse(status_code=500, content={"detail": str(e)})

    return FileResponse(path=png_path, media_type="image/png")


# ---------------------------------------------------------------------------
# Background threading event execution workers
# ---------------------------------------------------------------------------


async def _run_job(job_id: str, file_bytes: bytes, params: WebParams) -> None:
    """Execute single parsing inside isolation threads safely recording outputs."""
    _jobs[job_id]["status"] = "running"
    loop = asyncio.get_event_loop()
    try:
        result = await loop.run_in_executor(
            None,
            run_single,
            file_bytes,
            params.sample_name,
            params.dpi,
            params.color,
            params.group,
        )
        _jobs[job_id].update(
            {
                "status": "done",
                "path": result.zip_path,
                "png_path": result.png_path,
                "tsv_path": result.tsv_path,
                "samples": result.samples,
                "pathways": result.pathways,
            }
        )
    except Exception as e:
        _jobs[job_id].update({"status": "error", "message": str(e)})


async def _run_job_multi(
    job_id: str, named_files: list[tuple[str, bytes]], params: WebParams
) -> None:
    """Execute collection serialization loops mapping structural tabular datasets."""
    _jobs[job_id]["status"] = "running"
    loop = asyncio.get_event_loop()
    try:
        result = await loop.run_in_executor(
            None,
            run_multi,
            named_files,
            params.dpi,
            params.color,
            params.group,
        )
        _jobs[job_id].update(
            {
                "status": "done",
                "path": result.zip_path,
                "png_path": result.png_path,
                "tsv_path": result.tsv_path,
                "samples": result.samples,
                "pathways": result.pathways,
            }
        )
    except Exception as e:
        _jobs[job_id].update({"status": "error", "message": str(e)})


# ---------------------------------------------------------------------------
# Blocking visualization core loop sub-engine
# ---------------------------------------------------------------------------


def _blocking_viz(
    tsv_path: str,
    plot_type: str,
    figsize: tuple[float, float] | None,
    dpi: int,
    heatmap_color: str,
    heatmap_group: bool,
    heatmap_sample_name: str,
    heatmap_dpi: int,
    title: str,
    title_fontsize: float,
    title_color: str,
    title_weight: FontWeight,
    title_style: FontStyle,
    background_color: str,
    xlabel: str,
    xlabel_fontsize: float,
    xlabel_color: str,
    xlabel_weight: FontWeight,
    xlabel_style: FontStyle,
    ylabel: str,
    ylabel_fontsize: float,
    ylabel_color: str,
    ylabel_weight: FontWeight,
    ylabel_style: FontStyle,
    xticks_fontsize: float,
    xticks_color: str,
    xticks_weight: FontWeight,
    xticks_style: FontStyle,
    xticks_rotation: float,
    xticks_ha: HorizontalAlignment,
    yticks_fontsize: float,
    yticks_color: str,
    yticks_weight: FontWeight,
    yticks_style: FontStyle,
    grid: bool,
    grid_linestyle: str,
    grid_alpha: float,
    cmap: str,
    cmap_range_min: int,
    cmap_range_max: int,
    sort_order: SortOrder,
    box_color: str,
    showfliers: bool,
    grid_color: str,
    grid_linewidth: float,
    threshold: float,
    node_size: float,
    node_color: str,
    node_edgecolors: str,
    node_linewidths: float,
    label_fontsize: float,
    label_color: str,
    label_weight: FontWeight,
    edge_cmap_name: str,
    cbar_size: float,
    pathways_selected: list[str],
    colors_selected: list[str],
    sample_order: list[str],
    fill_alpha: float,
    line_width: float,
    line_style: LineStyle,
    label_background: str | None,
    label_edgecolor: str | None,
    label_pad: float,
    show_legend: bool,
    legend_loc: str,
    bar_width: float,
    edgecolor: str,
    edge_linewidth: float,
    stream_fill_alpha: float,
    legend_fontsize: float,
    legend_bbox: tuple[float, float],
) -> str:
    """Thread-isolated blocking core drawing dispatcher rendering static PNGs to temporary storage."""
    import matplotlib.pyplot as plt

    from .cheatmaps import (
        grouped_heatmap,
        grouped_heatmap_multi,
        simple_heatmap,
        simple_heatmap_multi,
    )
    from .kgnplot import barplot as bp
    from .kgnplot import corrnet as cn
    from .kgnplot import radarplot as rp
    from .kgnplot import stackedbar as sb
    from .kgnplot import streamgraph as sg

    df = pd.read_csv(tsv_path, sep="\t")

    if sample_order:
        ordered_cols = ["Function"] + [s for s in sample_order if s in df.columns]
        remaining = [c for c in df.columns if c not in ordered_cols]
        df = df[ordered_cols + remaining]

    title_val = title if title else None
    cmap_val = cmap if cmap else None
    out_path = _secure_temp_path(".png")

    if plot_type == "barplot":
        fs = figsize or (8.0, 12.0)
        plot = bp.barplot(
            df,
            figsize=fs,
            cmap=cmap_val or "Greens",
            cmap_range=(cmap_range_min, cmap_range_max),
            title=title_val,
            title_fontsize=title_fontsize,
            title_color=title_color,
            title_weight=title_weight,
            title_style=title_style,
            xlabel=xlabel or "Pathway completeness",
            xlabel_fontsize=xlabel_fontsize,
            xlabel_color=xlabel_color,
            xlabel_weight=xlabel_weight,
            xlabel_style=xlabel_style,
            ylabel=ylabel or "Pathway",
            ylabel_fontsize=ylabel_fontsize,
            ylabel_color=ylabel_color,
            ylabel_weight=ylabel_weight,
            ylabel_style=ylabel_style,
            xticks_fontsize=xticks_fontsize,
            xticks_color=xticks_color,
            xticks_weight=xticks_weight,
            xticks_style=xticks_style,
            yticks_fontsize=yticks_fontsize,
            yticks_color=yticks_color,
            yticks_weight=yticks_weight,
            yticks_style=yticks_style,
            grid=grid,
            grid_linestyle=grid_linestyle,
            grid_alpha=grid_alpha,
            background_color=background_color,
            sort_order=sort_order,
        )

    elif plot_type == "corrnet":
        import re

        import matplotlib.cm as mcm

        fs = figsize or (12.0, 6.0)
        cmap_key = edge_cmap_name or "coolwarm"
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", cmap_key) and hasattr(mcm, cmap_key):
            edge_cmap_obj = getattr(mcm, cmap_key)
        else:
            edge_cmap_obj = colormaps["coolwarm"]
        plot = cn.correlation_network(
            df,
            figsize=fs,
            threshold=threshold,
            node_size=node_size,
            node_color=node_color,
            node_edgecolors=node_edgecolors,
            node_linewidths=node_linewidths,
            label_fontsize=label_fontsize,
            label_color=label_color,
            label_weight=label_weight,
            edge_cmap=edge_cmap_obj,
            cbar_size=cbar_size,
            title=title_val,
            title_fontsize=title_fontsize,
            title_color=title_color,
            title_weight=title_weight,
            title_style=title_style,
            background_color=background_color,
        )

    elif plot_type == "radarplot":
        fs = figsize or (8.0, 8.0)
        if not pathways_selected:
            raise ValueError("Please select 1–4 pathways for the radar plot.")
        plot = rp.radarplot(
            df,
            pathways=pathways_selected,
            figsize=fs,
            colors=colors_selected or None,
            sample_order=sample_order or None,
            title=title_val,
            title_fontsize=title_fontsize,
            title_color=title_color,
            title_weight=title_weight,
            title_style=title_style,
            label_fontsize=label_fontsize,
            label_color=label_color,
            label_weight=label_weight,
            label_style=xticks_style,
            label_background=label_background,
            label_edgecolor=label_edgecolor,
            label_pad=label_pad,
            ytick_fontsize=yticks_fontsize,
            ytick_color=yticks_color,
            ytick_weight=yticks_weight,
            fill_alpha=fill_alpha,
            line_width=line_width,
            line_style=line_style,
            background_color=background_color,
            legend_loc=legend_loc,
            legend_bbox=legend_bbox,
            show_legend=show_legend,
        )

    elif plot_type == "stackedbar":
        fs = figsize or (14.0, 7.0)
        plot = sb.stacked_barplot(
            df,
            figsize=fs,
            cmap=cmap_val or "tab20",
            bar_width=bar_width,
            edgecolor=edgecolor,
            edge_linewidth=edge_linewidth,
            title=title_val,
            title_fontsize=title_fontsize,
            title_color=title_color,
            title_weight=title_weight,
            title_style=title_style,
            xlabel=xlabel or "Samples",
            xlabel_fontsize=xlabel_fontsize,
            xlabel_color=xlabel_color,
            xlabel_weight=xlabel_weight,
            xlabel_style=xlabel_style,
            ylabel=ylabel or "Total Completeness",
            ylabel_fontsize=ylabel_fontsize,
            ylabel_color=ylabel_color,
            ylabel_weight=ylabel_weight,
            ylabel_style=ylabel_style,
            xticks_rotation=xticks_rotation,
            xticks_ha=xticks_ha,
            xticks_fontsize=xticks_fontsize,
            xticks_color=xticks_color,
            xticks_weight=xticks_weight,
            xticks_style=xticks_style,
            background_color=background_color,
            grid=grid,
            grid_linestyle=grid_linestyle,
            grid_alpha=grid_alpha,
            legend_fontsize=legend_fontsize,
            legend_loc=legend_loc,
            legend_bbox=legend_bbox,
            show_legend=show_legend,
        )

    elif plot_type == "streamgraph":
        fs = figsize or (14.0, 7.0)
        plot = sg.streamgraph(
            df,
            figsize=fs,
            cmap=cmap_val or "tab20",
            bar_width=bar_width,
            fill_alpha=stream_fill_alpha,
            edgecolor=edgecolor or None,
            edge_linewidth=edge_linewidth,
            title=title_val,
            title_fontsize=title_fontsize,
            title_color=title_color,
            title_weight=title_weight,
            title_style=title_style,
            xlabel=xlabel or "Samples",
            xlabel_fontsize=xlabel_fontsize,
            xlabel_color=xlabel_color,
            xlabel_weight=xlabel_weight,
            xlabel_style=xlabel_style,
            ylabel=ylabel or "Total Completeness",
            ylabel_fontsize=ylabel_fontsize,
            ylabel_color=ylabel_color,
            ylabel_weight=ylabel_weight,
            ylabel_style=ylabel_style,
            xticks_rotation=xticks_rotation,
            xticks_ha=xticks_ha,
            xticks_fontsize=xticks_fontsize,
            xticks_color=xticks_color,
            xticks_weight=xticks_weight,
            xticks_style=xticks_style,
            background_color=background_color,
            grid=grid,
            grid_linestyle=grid_linestyle,
            grid_alpha=grid_alpha,
            legend_fontsize=legend_fontsize,
            legend_loc=legend_loc,
            legend_bbox=legend_bbox,
            show_legend=show_legend,
        )

    elif plot_type == "heatmap":
        is_multi = "Function" in df.columns
        try:
            with tempfile.TemporaryDirectory() as hm_tmp:
                hm_tmp_path = Path(hm_tmp)
                if is_multi:
                    if heatmap_group:
                        grouped_heatmap_multi.generate_grouped_heatmap_multi(
                            df, str(hm_tmp_path), heatmap_dpi, heatmap_color
                        )
                    else:
                        simple_heatmap_multi.generate_heatmap_multi(
                            df, str(hm_tmp_path), heatmap_dpi, heatmap_color
                        )
                else:
                    sample_name = (
                        heatmap_sample_name if heatmap_sample_name else "SAMPLE"
                    )
                    if heatmap_group:
                        grouped_heatmap.generate_grouped_heatmap(
                            tsv_path,
                            str(hm_tmp_path),
                            heatmap_dpi,
                            heatmap_color,
                            sample_name,
                        )
                    else:
                        simple_heatmap.generate_heatmap(
                            tsv_path,
                            str(hm_tmp_path),
                            heatmap_dpi,
                            heatmap_color,
                            sample_name,
                        )
                png_files = list(hm_tmp_path.rglob("*.png"))
                if not png_files:
                    raise FileNotFoundError("Heatmap generation produced no PNG.")
                shutil.copy2(png_files[0], out_path)
        finally:
            plt.close("all")
        return str(out_path)

    else:
        raise ValueError(f"Unknown plot_type: {plot_type!r}")

    plot.savefig(str(out_path), dpi=dpi)
    plt.close("all")
    return str(out_path)
