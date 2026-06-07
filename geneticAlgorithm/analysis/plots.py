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
    text = (
        f"Generations recorded: {len(generations)}  |  "
        f"Best fitness: {max(best_fitness):.4f}  |  "
        f"Final average fitness: {average_fitness[-1]:.4f}"
    )
    if len(best_fitness) >= 2:
        gain = best_fitness[-1] - best_fitness[0]
        text += f"  |  Best fitness gain: {gain:+.4f}"
    return text


def _selection_label(config):
    if config.get("roulette_selection"):
        return "Roulette wheel"
    return f"Tournament (size {config.get('tournament_size')})"


def run_config_table_rows(config, initial_population_size=None):
    if not config:
        return []

    rows = [
        ("Planned generations", str(config.get("generations", "n/a"))),
        ("Selection", _selection_label(config)),
        ("Elitism", str(config.get("elitism_size", "n/a"))),
        ("Mutation probability", str(config.get("mutation_probability", "n/a"))),
        ("Crossover mode", str(config.get("crossover_mode", "n/a"))),
        ("Mutation mode", str(config.get("mutation_mode", "n/a"))),
    ]
    if initial_population_size is not None:
        rows.insert(4, ("Initial population", str(initial_population_size)))
    return rows


def run_config_summary_text(config, initial_population_size=None):
    rows = run_config_table_rows(config, initial_population_size=initial_population_size)
    if not rows:
        return "Run configuration not available."
    return "  |  ".join(f"{label}: {value}" for label, value in rows)


def _format_rate(summary, key):
    value = summary.get(key) if summary else None
    if value is None:
        return "n/a"
    return f"{value:.1f}%"


def operators_table_rows(operators_summary):
    if not operators_summary:
        return []

    crossover = operators_summary.get("crossover", {})
    mutation = operators_summary.get("mutation", {})
    crossover_stats = crossover.get("stats", {})
    mutation_stats = mutation.get("stats", {})
    crossover_summary = crossover.get("summary", {})
    mutation_summary = mutation.get("summary", {})

    rows = [
        ("__section__", "Crossover"),
        ("Mode", str(crossover.get("mode", "n/a"))),
        (
            "Successes",
            f"{crossover_stats.get('successes', 0)} / {crossover_stats.get('attempts', 0)}",
        ),
        ("Success rate", _format_rate(crossover_summary, "success_rate_pct")),
        ("Failure rate", _format_rate(crossover_summary, "failure_rate_pct")),
        ("__section__", "Mutation"),
        ("Mode", str(mutation.get("mode", "n/a"))),
        (
            "Successes",
            f"{mutation_stats.get('successes', 0)} / {mutation_stats.get('invocations', 0)}",
        ),
        ("Success rate", _format_rate(mutation_summary, "success_rate_pct")),
        ("Fallback rate", _format_rate(mutation_summary, "fallback_rate_pct")),
        (
            "Avg attempts on success",
            str(mutation_summary.get("avg_attempts_on_success", "n/a")),
        ),
    ]

    by_kind = mutation_stats.get("by_kind") or {}
    if by_kind:
        rows.append(("__section__", "Mutation by kind"))
        for kind, count in sorted(by_kind.items()):
            rows.append((kind, str(count)))

    return rows


def operators_summary_text(operators_summary):
    rows = operators_table_rows(operators_summary)
    if not rows:
        return (
            "Operator statistics not available yet. "
            "Complete at least one full GA run (Launch or Generate)."
        )

    lines = []
    for label, value in rows:
        if label == "__section__":
            lines.append(value)
        else:
            lines.append(f"{label}: {value}")
    return "\n".join(lines)


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


def best_molecule_summary_text(best_molecule):
    if not best_molecule:
        return "Best molecule not available."

    lines = []
    name = (best_molecule.get("description") or "").strip()
    if name:
        lines.append(f"Name: {name}")
    lines.append(f"Fitness: {best_molecule.get('fitness', 0.0):.4f}")
    generation = best_molecule.get("generation")
    if generation is not None:
        lines.append(f"Found in generation: {generation}")
    smiles = best_molecule.get("smiles") or ""
    if smiles:
        lines.append(f"SMILES: {smiles}")
    return "\n".join(lines)


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
