"""Shared plotting helpers for evolution run analysis (app + notebooks)."""

import statistics

from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator


def plot_fitness_over_generations(ax, generations, best_fitness, average_fitness):
    ax.plot(
        generations,
        best_fitness,
        marker="o",
        markersize=3,
        linewidth=1.4,
        color="#2e7d32",
        label="Best fitness",
    )
    ax.plot(
        generations,
        average_fitness,
        marker="o",
        markersize=3,
        linewidth=1.4,
        color="#1565c0",
        label="Average fitness",
    )
    ax.set_xlabel("Generation", fontsize=8)
    ax.set_ylabel("Fitness", fontsize=8)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(labelsize=7)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(loc="best", fontsize=7)


def create_fitness_figure(generations, best_fitness, average_fitness, figsize=(3.0, 1.7)):
    figure = Figure(figsize=figsize)
    axes = figure.add_subplot(111)
    plot_fitness_over_generations(axes, generations, best_fitness, average_fitness)
    figure.tight_layout(pad=1.0)
    return figure


def fitness_summary_text(generations, best_fitness, average_fitness):
    if not generations:
        return ""
    return (
        f"Generations recorded: {len(generations)}  |  "
        f"Best fitness: {max(best_fitness):.4f}  |  "
        f"Final average fitness: {average_fitness[-1]:.4f}"
    )


def plot_tanimoto_histogram(ax, similarities, generation=None, bins=20):
    if not similarities:
        ax.text(
            0.5,
            0.5,
            "Not enough valid molecules\nfor pairwise comparison",
            ha="center",
            va="center",
            fontsize=8,
            transform=ax.transAxes,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        return

    ax.hist(
        similarities,
        bins=bins,
        range=(0.0, 1.0),
        color="#5c6bc0",
        edgecolor="white",
        linewidth=0.5,
    )
    if generation is not None:
        ax.set_title(f"Generation {generation}", fontsize=9)
    ax.set_xlabel("Tanimoto similarity", fontsize=8)
    ax.set_ylabel("Pair count", fontsize=8)
    ax.set_xlim(0.0, 1.0)
    ax.tick_params(labelsize=7)
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)


def create_tanimoto_histogram_figure(similarities, generation=None, figsize=(3.0, 1.7)):
    figure = Figure(figsize=figsize)
    axes = figure.add_subplot(111)
    plot_tanimoto_histogram(axes, similarities, generation=generation)
    figure.tight_layout(pad=1.0)
    return figure


def tanimoto_summary_text(diversity_summary, generation=None, similarities=None):
    if not diversity_summary or diversity_summary.get("pair_count", 0) == 0:
        return "No pairwise Tanimoto similarities available."

    std = diversity_summary.get("std")
    if std is None and similarities:
        std = round(statistics.pstdev(similarities), 4)

    gen_prefix = f"Generation {generation}: " if generation is not None else ""
    std_text = f"std: {std:.4f}  |  " if std is not None else ""
    return (
        f"{gen_prefix}"
        f"{diversity_summary['pair_count']} pairs  |  "
        f"mean: {diversity_summary['mean']:.4f}  |  "
        f"{std_text}"
        f"min: {diversity_summary['min']:.4f}  |  "
        f"max: {diversity_summary['max']:.4f}"
    )
