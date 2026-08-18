#!/usr/bin/env python3
"""Plot mean best/average fitness curves for population-size experiments."""

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
POPULATION_SIZES = [15, 30, 50, 100]
OUTPUT_DIR = Path("results/plots")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RUN_DIR_PATTERN = re.compile(r"^gen(\d+)_run(\d+)$", re.IGNORECASE)


def discover_runs(runs_dir: Path) -> dict[int, list[Path]]:
    """Group generations.json files by population size."""
    grouped: dict[int, list[Path]] = defaultdict(list)

    for run_dir in sorted(runs_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        match = RUN_DIR_PATTERN.match(run_dir.name)
        if not match:
            continue
        pop_size = int(match.group(1))
        json_path = run_dir / "generations.json"
        if json_path.is_file():
            grouped[pop_size].append(json_path)

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
    curves: dict[int, tuple[list[int], list[float]]],
    output_stem: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    colors = {
        15: "#2e7d32",
        30: "#1565c0",
        50: "#ef6c00",
        100: "#6a1b9a",
    }

    for pop_size in POPULATION_SIZES:
        if pop_size not in curves:
            print(f"[WARN] Nema podataka za populaciju {pop_size}")
            continue

        xs, ys = curves[pop_size]
        ax.plot(
            xs,
            ys,
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=f"Популација {pop_size}",
            color=colors.get(pop_size),
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

    best_curves: dict[int, tuple[list[int], list[float]]] = {}
    avg_curves: dict[int, tuple[list[int], list[float]]] = {}

    for pop_size in POPULATION_SIZES:
        paths = grouped.get(pop_size, [])
        if not paths:
            continue

        print(f"\Популација {pop_size}: {len(paths)} извршавања")
        for p in paths:
            print(f"  - {p.parent.name}")

        best_runs = []
        avg_runs = []

        for path in paths:
            gens, best, avg = extract_series(path)
            best_runs.append((gens, best))
            avg_runs.append((gens, avg))

        best_curves[pop_size] = mean_series(best_runs)
        avg_curves[pop_size] = mean_series(avg_runs)

    plot_metric(
        metric_name="Најбољи фитнес",
        ylabel="Најбољи фитнес",
        curves=best_curves,
        output_stem="best_fitness_by_population",
    )
    plot_metric(
        metric_name="Просечан фитнес",
        ylabel="Просечан фитнес",
        curves=avg_curves,
        output_stem="average_fitness_by_population",
    )


if __name__ == "__main__":
    main()