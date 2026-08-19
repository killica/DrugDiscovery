#!/usr/bin/env python3
"""Plot mean best/average fitness curves for genetic-operator experiments."""

from __future__ import annotations

import re
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from analysis.io import load_generations_export

# --- configuration ---
RUNS_DIR = Path("results/runs")
OUTPUT_DIR = Path("results/plots")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OPERATOR_CONFIGS = ["SMILES", "SELFIES", "GRAPH", "BRICS"]
RUN_DIR_PATTERN = re.compile(r"^op_(smiles|selfies|graph|brics)_run(\d+)$", re.IGNORECASE)

LABELS = {
    "SMILES": "SMILES",
    "SELFIES": "SELFIES",
    "GRAPH": "Граф",
    "BRICS": "BRICS",
}

COLORS = {
    "SMILES": "#1565c0",
    "SELFIES": "#2e7d32",
    "GRAPH": "#ef6c00",
    "BRICS": "#6a1b9a",
}


def discover_runs(runs_dir: Path) -> dict[str, list[Path]]:
    """Group generations.json files by operator configuration."""
    grouped: dict[str, list[Path]] = defaultdict(list)

    for run_dir in sorted(runs_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        match = RUN_DIR_PATTERN.match(run_dir.name)
        if not match:
            continue
        config = match.group(1).upper()
        json_path = run_dir / "generations.json"
        if json_path.is_file():
            grouped[config].append(json_path)

    for config in grouped:
        grouped[config].sort(key=lambda p: p.parent.name)

    return grouped


def extract_series(path: Path) -> tuple[list[int], list[float], list[float]]:
    stats = load_generations_export(str(path))
    generations = stats.generation_indices
    return generations, stats.best_fitness, stats.average_fitness


def mean_series(all_series: list[tuple[list[int], list[float]]]) -> tuple[list[int], list[float]]:
    """Average y-values across runs for each generation index."""
    if not all_series:
        return [], []

    n_generations = min(len(gens) for gens, _ in all_series)
    generations = all_series[0][0][:n_generations]

    means = []
    for i in range(n_generations):
        values = [ys[i] for _, ys in all_series if i < len(ys)]
        means.append(statistics.mean(values))

    return generations, means


def plot_metric(
    metric_name: str,
    ylabel: str,
    curves: dict[str, tuple[list[int], list[float]]],
    output_stem: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    for config in OPERATOR_CONFIGS:
        if config not in curves:
            print(f"[WARN] Nema podataka za konfiguraciju {config}")
            continue

        xs, ys = curves[config]
        ax.plot(
            xs,
            ys,
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=LABELS[config],
            color=COLORS.get(config),
        )

    ax.set_xlabel("Генерација")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{metric_name} кроз генерације (средња вредност)")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(loc="best")
    fig.tight_layout()

    png_path = OUTPUT_DIR / f"{output_stem}.png"
    pdf_path = OUTPUT_DIR / f"{output_stem}.pdf"
    fig.savefig(png_path, dpi=200)
    fig.savefig(pdf_path)
    plt.close(fig)

    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


def main() -> None:
    grouped = discover_runs(RUNS_DIR)

    best_curves: dict[str, tuple[list[int], list[float]]] = {}
    avg_curves: dict[str, tuple[list[int], list[float]]] = {}

    for config in OPERATOR_CONFIGS:
        paths = grouped.get(config, [])
        if not paths:
            continue

        print(f"\n{LABELS[config]}: {len(paths)} izvršavanja")
        for p in paths:
            print(f"  - {p.parent.name}")

        best_runs = []
        avg_runs = []

        for path in paths:
            gens, best, avg = extract_series(path)
            best_runs.append((gens, best))
            avg_runs.append((gens, avg))

        best_curves[config] = mean_series(best_runs)
        avg_curves[config] = mean_series(avg_runs)

    plot_metric(
        metric_name="Најбољи фитнес",
        ylabel="Најбољи фитнес",
        curves=best_curves,
        output_stem="best_fitness_by_operator",
    )
    plot_metric(
        metric_name="Просечни фитнес",
        ylabel="Просечни фитнес",
        curves=avg_curves,
        output_stem="average_fitness_by_operator",
    )


if __name__ == "__main__":
    main()