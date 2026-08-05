"""Pydantic schemas for data validation in KEGGaNOG.

This module defines the structural contracts used across both the
command-line interface (CLI) and the FastAPI web service, ensuring
consistent input validation and type safety for analysis parameters.
"""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

# ---------------------------------------------------------------------------
# Allowed values
# ---------------------------------------------------------------------------

ValidColor = Literal["Blues", "Greens", "Reds", "Purples", "Greys", "Oranges"]


# ---------------------------------------------------------------------------
# Base
# ---------------------------------------------------------------------------


class BaseParams(BaseModel):
    """Shared parameters for both CLI and web form."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    dpi: int = Field(default=300, ge=72, le=600)
    color: ValidColor = Field(default="Blues")
    sample_name: str = Field(default="SAMPLE", max_length=64)
    group: bool = Field(default=False)

    @field_validator("sample_name")
    @classmethod
    def no_path_characters(cls, v: str) -> str:
        """Reject sample names with filesystem-unsafe characters."""
        forbidden = set(r'\/:*?"<>|')
        bad_chars = forbidden & set(v)
        if bad_chars:
            raise ValueError(f"sample_name contains invalid characters: {bad_chars}")
        return v


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class CLIParams(BaseParams):
    """Validated parameters for the KEGGaNOG command-line interface."""

    input_path: str = Field(..., description="Path to input file/list.")
    output_dir: str = Field(..., description="Output folder.")
    multi: bool = Field(default=False)
    overwrite: bool = Field(default=False)


class WebParams(BaseParams):
    """Validated parameters for the KEGGaNOG web form."""


class JobStatus(BaseModel):
    """Status of a background analysis job."""

    model_config = ConfigDict(frozen=True)

    job_id: str
    status: Literal["pending", "running", "done", "error"]
    message: str = Field(default="")
