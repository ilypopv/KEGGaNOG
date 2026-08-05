"""Main command-line interface entry point for KEGGaNOG.

This module orchestrates the entire KEGGaNOG suite, exposing a unified interface
for linking eggNOG-mapper functional annotations with KEGG-Decoder metabolic
pathway reconstructions, supporting both automated workflows and a local web UI.
"""

import logging
import shutil
import sys
from importlib.metadata import version as _metadata_version
from pathlib import Path

import typer
from pydantic import ValidationError

from . import kegganog_multi
from .processing.pipeline import run_single
from .schemas import CLIParams, ValidColor

# Initialize module-level isolated logger
_log: logging.Logger = logging.getLogger(__name__)

# Primary Typer application instantiation
app: typer.Typer = typer.Typer(
    name="kegganog",
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)

# ===========================================================================
# Academic Disclosures & Citation
# ===========================================================================

CITATION_MESSAGE: str = """
Thank you for using KEGGaNOG! If you use it in your research, please cite:

    Popov, I.V., Chikindas, M.L., Venema, K., Ermakov, A.M. and Popov, I.V., 2025.
    KEGGaNOG: A Lightweight Tool for KEGG Module Profiling From Orthology-Based Annotations.
    Molecular Nutrition & Food Research, p.e70269.
    doi.org/10.1002/mnfr.70269
"""


def print_citation() -> None:
    """Print the official academic citation message to standard output."""
    print(CITATION_MESSAGE)


# ===========================================================================
# Internal Helpers & Operational Callbacks
# ===========================================================================


def _version_callback(value: bool) -> None:
    """Eager callback executing the version flag logic.

    Args:
        value: Boolean trigger provided by the Typer parameter evaluation.

    Raises:
        typer.Exit: Gracefully terminates runtime execution upon displaying version.
    """
    if value:
        print(f"KEGGaNOG {_metadata_version('kegganog')}")
        raise typer.Exit()


def _validate_params(
    input_path: str | None,
    output_dir: str | None,
    multi: bool,
    overwrite: bool,
    dpi: int,
    color: ValidColor,
    sample_name: str,
    group: bool,
) -> CLIParams:
    """Validate incoming CLI parameters against the Pydantic data engine schema.

    Args:
        input_path: Source path targeting eggNOG-mapper output tables.
        output_dir: Destination folder path for generated tables and figures.
        multi: Flag switching the engine matrix profile run state to multi-sample.
        overwrite: Overwrite protection bypass toggle.
        dpi: Target pixel density for visualization layers.
        color: Valid Seaborn color map palette identifier literal.
        sample_name: Character string tag representing the sample.
        group: Flag indicating categorical block clustering for pathways.

    Returns:
        CLIParams: A validated data container instance.

    Raises:
        typer.Exit: Terminated with exit code 1 if parsing constraints fail.
    """
    # Guarantee to static analyzers that optional paths have been checked and are not None
    assert input_path is not None
    assert output_dir is not None

    try:
        return CLIParams(
            input_path=input_path,
            output_dir=output_dir,
            multi=multi,
            overwrite=overwrite,
            dpi=dpi,
            color=color,
            sample_name=sample_name,
            group=group,
        )
    except ValidationError as e:
        print("KEGGaNOG: invalid argument value(s):", file=sys.stderr)
        for error in e.errors():
            field = error["loc"][0]
            message = error["msg"]
            print(f"  --{field}: {message}", file=sys.stderr)
        raise typer.Exit(code=1)


def _prepare_output_dir(params: CLIParams) -> None:
    """Initialize destination directory structures and wipe old records if forced.

    Args:
        params: Validated runtime parameter container context.

    Raises:
        FileExistsError: Overwrite protection checkpoint if target space is occupied.
    """
    output_dir: Path = Path(params.output_dir)

    if output_dir.exists():
        if not params.overwrite:
            raise FileExistsError(
                f"Output directory '{params.output_dir}' already exists. "
                "Use --overwrite to bypass overwrite protection."
            )
        shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True)
    (output_dir / "temp_files").mkdir()


def _run_pipeline(params: CLIParams) -> None:
    """Execute the single-sample metabolic pathway extraction pipeline.

    Args:
        params: Validated configuration schema object instance.
    """
    input_path: Path = Path(params.input_path)
    file_bytes: bytes = input_path.read_bytes()

    run_single(
        file_bytes=file_bytes,
        sample_name=params.sample_name,
        dpi=params.dpi,
        color=params.color,
        group=params.group,
        output_dir=params.output_dir,
    )


# ===========================================================================
# Primary Core CLI Command
# ===========================================================================


@app.command(
    help="KEGGaNOG: Link eggNOG-mapper and KEGG-Decoder for pathway visualization.",
    no_args_is_help=True,
)
def main(
    input_path: str | None = typer.Option(
        None,
        "-i",
        "--input",
        help="Path to eggNOG-mapper annotation file.",
    ),
    output_dir: str | None = typer.Option(
        None,
        "-o",
        "--output",
        help="Output folder to save results.",
    ),
    multi: bool = typer.Option(
        False,
        "-M",
        "--multi",
        help="Run KEGGaNOG in multi-sample cohort profile mode.",
    ),
    overwrite: bool = typer.Option(
        False,
        "-overwrite",
        "--overwrite",
        help="Overwrite the output directory if it already exists.",
    ),
    dpi: int = typer.Option(
        300,
        "-dpi",
        "--dpi",
        help="DPI resolution mapping index for the output image visualization.",
    ),
    color: ValidColor = typer.Option(
        "Blues",
        "-c",
        "--color",
        help="Target seaborn color map palette matrix rule.",
    ),
    sample_name: str = typer.Option(
        "SAMPLE",
        "-n",
        "--name",
        help="Sample identity text string for axis labeling.",
    ),
    group: bool = typer.Option(
        False,
        "-g",
        "--group",
        help="Group the pathway matrix heatmap rows based on predefined functional categories.",
    ),
    web: bool = typer.Option(
        False,
        "--web",
        help="Launch the interactive local web user interface dashboard at http://localhost:8000.",
    ),
    version: bool | None = typer.Option(
        None,
        "-V",
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    print("KEGGaNOG by Ilia V. Popov")

    # Intercept execution to delegate control to the Web Engine UI layer if requested
    if web:
        from .web import launch

        launch()
        return

    # Ensure CLI mandatory parameters are defined
    if input_path is None or output_dir is None:
        print(
            "Error: Missing mandatory options '--input' and '--output' in CLI production mode.",
            file=sys.stderr,
        )
        raise typer.Exit(code=1)

    # Parse parameter configurations through structural typing schema validation
    params: CLIParams = _validate_params(
        input_path=input_path,
        output_dir=output_dir,
        multi=multi,
        overwrite=overwrite,
        dpi=dpi,
        color=color,
        sample_name=sample_name,
        group=group,
    )

    # Setup local IO environment targets securely
    try:
        _prepare_output_dir(params)
    except FileExistsError as e:
        print(f"Error: {e}", file=sys.stderr)
        raise typer.Exit(code=1)

    # Route control context to specialized Single-Sample or Multi-Cohort runner layers
    if params.multi:
        assert params.input_path is not None
        assert params.output_dir is not None

        kegganog_multi.MultiSampleRunner(
            input_path=params.input_path,
            output_dir=params.output_dir,
            dpi=params.dpi,
            color=params.color,
            group=params.group,
        ).run()
    else:
        _run_pipeline(params)

    # Print citation contracts to the user
    print_citation()


def entry_point() -> None:
    """Consolidated system application execution manager.

    Intercepts runtime execution interrupts and secure context teardowns cleanly
    to maintain system trace sanity.
    """
    try:
        app()
    except KeyboardInterrupt:
        print("\nExecution runtime manually interrupted by user.", file=sys.stderr)
        print_citation()
        sys.exit(1)


if __name__ == "__main__":
    entry_point()
