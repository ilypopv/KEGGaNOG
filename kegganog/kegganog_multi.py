"""Orchestration engine for multi-sample functional profile processing.

This module coordinates multi-sample data ingestion from file lists or single
targets, parses multiple eggNOG-mapper output records, and coordinates aggregate
matrix reconstructions and cohort pathway visualization layouts.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from .processing.pipeline import PipelineResult, run_multi
from .schemas import ValidColor

# Initialize module-level isolated logger
_log: logging.Logger = logging.getLogger(__name__)


@dataclass
class MultiSampleRunner:
    """Orchestrates the multi-sample cohort analytical execution pipeline.

    Resolves collections of text-based sample paths from individual entries or
    manifest list files, coordinates disk-read transactions, and handles aggregate
    functional annotation evaluation.

    Attributes:
        input_path: Path targeting a single annotation document or a text file manifest
            containing one valid file path pointer per line.
        output_dir: Destination folder directory where files and summaries are saved.
        dpi: Pixel density mapping metric allocated to downstream visual charts.
        color: Target valid seaborn color map palette matrix rule literal.
        group: Flag indicating categorical grouping block clustering for pathways.
    """

    input_path: str
    output_dir: str
    dpi: int = 300
    color: ValidColor = "Blues"
    group: bool = False

    def run(self) -> PipelineResult:
        """Execute the complete multi-sample comparative metabolic matrix pipeline.

        Returns:
            PipelineResult: NamedTuple container hosting multi-sample summary matrix structures.
        """
        named_files: list[tuple[str, bytes]] = self._load_files()

        return run_multi(
            named_files=named_files,
            dpi=self.dpi,
            color=self.color,
            group=self.group,
            output_dir=self.output_dir,
        )

    def _collect_input_paths(self) -> list[str]:
        """Resolve individual annotation table paths from single targets or text manifests.

        Returns:
            list[str]: Array list pointing to raw string path elements.
        """
        # Parse multiline entry documents if a text manifest format is detected
        if self.input_path.endswith(".txt"):
            with open(self.input_path, encoding="utf-8") as fh:
                return [line.strip() for line in fh if line.strip()]

        # Fallback to single targeted extraction paths
        return [self.input_path]

    def _load_files(self) -> list[tuple[str, bytes]]:
        """Read target file segments into atomic bytes arrays from disk storage blocks.

        Skips missing file elements securely, appending descriptive diagnostic
        warnings to active log stream handlers.

        Returns:
            list[tuple[str, bytes]]: Sequence of resolved paired tuples (filename, content_bytes).
        """
        named_files: list[tuple[str, bytes]] = []

        for file_path in self._collect_input_paths():
            fp: Path = Path(file_path)

            # Atomic safety validation step protecting IO transactions
            if not fp.is_file():
                _log.warning(
                    "Input file %s does not exist on disk. Skipping.", file_path
                )
                continue

            named_files.append((fp.name, fp.read_bytes()))

        return named_files
