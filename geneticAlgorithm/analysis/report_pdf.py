"""Export evolution run statistics to a PDF report."""

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from rdkit import Chem
from rdkit.Chem import Draw

from analysis.plots import (
    plot_fitness_over_generations,
    plot_tanimoto_histogram,
    fitness_summary_text,
    run_config_table_rows,
    operators_table_rows,
    best_molecule_summary_text,
    tanimoto_summary_text,
)


def _table_cell_rows(rows):
    cells = []
    for label, value in rows:
        if label == "__section__":
            cells.append([str(value), ""])
        else:
            cells.append([str(label), str(value)])
    return cells


def _draw_parameter_table(ax, title, rows):
    ax.axis("off")
    ax.set_title(title, fontsize=10, fontweight="bold", pad=8)
    if not rows:
        ax.text(0.5, 0.5, "No data available.", ha="center", va="center", fontsize=9)
        return

    table = ax.table(
        cellText=_table_cell_rows(rows),
        colLabels=["Parameter", "Value"],
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.15)

    for row_index, (label, _value) in enumerate(rows):
        table_row = row_index + 1
        if label == "__section__":
            for col in (0, 1):
                cell = table[(table_row, col)]
                cell.set_facecolor("#e8f5e9")
                cell.get_text().set_fontweight("bold")
                cell.get_text().set_color("#1b5e20")


def _draw_best_molecule_panel(ax, best_molecule, image_size=250):
    ax.axis("off")
    ax.set_title("Best molecule so far", fontsize=10, fontweight="bold", pad=8)

    if not best_molecule or not best_molecule.get("smiles"):
        ax.text(0.5, 0.5, "Best molecule not available.", ha="center", va="center", fontsize=9)
        return

    mol = Chem.MolFromSmiles(best_molecule["smiles"])
    if mol is not None:
        mol_image = Draw.MolToImage(mol, size=(image_size, image_size))
        ax.imshow(mol_image, aspect="equal")
        ax.set_anchor("N")

    details = best_molecule_summary_text(best_molecule)
    ax.text(
        0.5,
        -0.08,
        details,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7.5,
        wrap=True,
    )


def export_evolution_report_pdf(stats, output_path, best_molecule=None):
    """Write the current evolution report to a PDF file."""
    if not stats.has_data():
        raise ValueError("No evolution data available to export.")

    if best_molecule is None:
        best_molecule = getattr(stats, "best_candidate", None)

    initial_population_size = None
    if stats.generations_data:
        initial_population_size = stats.generations_data[0].get("population_size")

    config_rows = run_config_table_rows(
        stats.config,
        initial_population_size=initial_population_size,
    )
    operator_rows = operators_table_rows(stats.operators_summary)

    final_entry = stats.final_generation_entry()
    diversity = final_entry.get("diversity") if final_entry else None
    final_generation = final_entry.get("generation") if final_entry else None

    with PdfPages(output_path) as pdf:
        summary_figure = Figure(figsize=(11.0, 8.5))
        grid = GridSpec(2, 3, figure=summary_figure, height_ratios=[1.35, 1.0], hspace=0.45, wspace=0.35)

        title_lines = ["Evolution report"]
        if stats.run_id:
            title_lines.append(f"Run ID: {stats.run_id}")
        if stats.started_at:
            title_lines.append(f"Started: {stats.started_at}")
        summary_figure.suptitle("\n".join(title_lines), fontsize=14, fontweight="bold", y=0.98)

        _draw_parameter_table(
            summary_figure.add_subplot(grid[0, 0]),
            "Run settings",
            config_rows,
        )
        _draw_best_molecule_panel(
            summary_figure.add_subplot(grid[0, 1]),
            best_molecule,
        )
        _draw_parameter_table(
            summary_figure.add_subplot(grid[0, 2]),
            "Genetic operators (full run)",
            operator_rows,
        )

        fitness_ax = summary_figure.add_subplot(grid[1, 0])
        plot_fitness_over_generations(
            fitness_ax,
            stats.generation_indices,
            stats.best_fitness,
            stats.average_fitness,
        )
        fitness_ax.set_title("Fitness over generations", fontsize=10, fontweight="bold")
        fitness_ax.text(
            0.5,
            -0.28,
            fitness_summary_text(
                stats.generation_indices,
                stats.best_fitness,
                stats.average_fitness,
            ),
            transform=fitness_ax.transAxes,
            ha="center",
            va="top",
            fontsize=7.5,
            wrap=True,
        )

        diversity_ax = summary_figure.add_subplot(grid[1, 1:])
        if diversity and diversity.get("summary", {}).get("pair_count", 0) > 0:
            plot_tanimoto_histogram(
                diversity_ax,
                diversity.get("pairwise_tanimoto", []),
                generation=final_generation,
            )
            diversity_ax.set_title(
                "Population diversity (final generation)",
                fontsize=10,
                fontweight="bold",
            )
            diversity_ax.text(
                0.5,
                -0.28,
                tanimoto_summary_text(
                    diversity.get("summary"),
                    generation=final_generation,
                    similarities=diversity.get("pairwise_tanimoto", []),
                ),
                transform=diversity_ax.transAxes,
                ha="center",
                va="top",
                fontsize=7.5,
                wrap=True,
            )
        else:
            diversity_ax.axis("off")
            diversity_ax.set_title(
                "Population diversity (final generation)",
                fontsize=10,
                fontweight="bold",
            )
            diversity_ax.text(
                0.5,
                0.5,
                "No pairwise Tanimoto similarities available.",
                ha="center",
                va="center",
                fontsize=9,
            )

        summary_figure.tight_layout(rect=(0, 0, 1, 0.94))
        pdf.savefig(summary_figure)
        summary_figure.clear()

    return output_path
