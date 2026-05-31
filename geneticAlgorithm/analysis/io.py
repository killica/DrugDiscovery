"""Load structured evolution run exports."""

import os

from evolutionRun import EvolutionStatistics, GENERATIONS_FILENAME


def load_generations_export(path):
    """Load a run from ``results/runs/<run_id>/generations.json``."""
    return EvolutionStatistics.from_export(path)


def load_latest_run_export(runs_dir="results/runs"):
    """Load the most recently modified ``generations.json`` under ``runs_dir``."""
    if not os.path.isdir(runs_dir):
        raise FileNotFoundError(f"No runs directory: {runs_dir}")

    candidates = []
    for run_id in os.listdir(runs_dir):
        path = os.path.join(runs_dir, run_id, GENERATIONS_FILENAME)
        if os.path.isfile(path):
            candidates.append((os.path.getmtime(path), path))

    if not candidates:
        raise FileNotFoundError(f"No {GENERATIONS_FILENAME} files in {runs_dir}")

    return EvolutionStatistics.from_export(max(candidates)[1])
