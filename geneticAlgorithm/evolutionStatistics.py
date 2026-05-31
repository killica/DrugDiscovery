from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QGroupBox,
    QGraphicsDropShadowEffect,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from analysis.plots import (
    plot_fitness_over_generations,
    plot_tanimoto_histogram,
    fitness_summary_text,
    tanimoto_summary_text,
)
from evolutionRun import EvolutionStatistics


STATS_CARD_QSS = """
    QGroupBox#statsCard {
        background-color: rgb(255, 255, 255);
        border: 2px solid #43a047;
        border-radius: 15px;
        padding: 14px;
        margin-top: 10px;
        font-size: 14px;
        font-weight: bold;
        color: #1b5e20;
    }
    QGroupBox#statsCard::title {
        subcontrol-origin: margin;
        left: 14px;
        padding: 0 8px;
    }
"""


def _apply_card_shadow(widget):
    shadow = QGraphicsDropShadowEffect(widget)
    shadow.setOffset(3, 5)
    shadow.setBlurRadius(20)
    shadow.setColor(QColor(0, 0, 0, 95))
    widget.setGraphicsEffect(shadow)


class _MatplotlibCardChart(FigureCanvas):
    CHART_WIDTH = 400
    CHART_HEIGHT = 250

    def __init__(self, parent=None, figsize=(3.0, 1.7)):
        figure = Figure(figsize=figsize)
        self._axes = figure.add_subplot(111)
        super().__init__(figure)
        self.setParent(parent)
        self.setFixedSize(self.CHART_WIDTH, self.CHART_HEIGHT)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        figure.tight_layout(pad=1.0)


class FitnessMatplotlibChart(_MatplotlibCardChart):
    def plot(self, generations, best_fitness, average_fitness):
        self._axes.clear()
        plot_fitness_over_generations(
            self._axes,
            generations,
            best_fitness,
            average_fitness,
        )
        self.figure.tight_layout(pad=1.0)
        self.draw()


class TanimotoHistogramChart(_MatplotlibCardChart):
    def plot(self, similarities, generation=None):
        self._axes.clear()
        plot_tanimoto_histogram(
            self._axes,
            similarities,
            generation=generation,
        )
        self.figure.tight_layout(pad=1.0)
        self.draw()


class EvolutionStatsChart(QWidget):
    """Scrollable report panel; add more chart sections here over time."""

    CARD_WIDTH = 440

    def __init__(self):
        super().__init__()
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setStyleSheet("QScrollArea { background: transparent; border: none; }")
        outer.addWidget(scroll)

        content = QWidget()
        self._content_layout = QVBoxLayout(content)
        self._content_layout.setContentsMargins(16, 16, 28, 28)
        self._content_layout.setSpacing(18)

        self.empty_label = QLabel(
            "No evolution data yet.\nRun the genetic algorithm to collect statistics."
        )
        self.empty_label.setAlignment(Qt.AlignCenter)
        self.empty_label.setStyleSheet("color: #616161; font-size: 14px;")
        self._content_layout.addWidget(self.empty_label)

        self.fitness_card, self.fitness_chart, self.fitness_summary_label = self._make_card(
            "Fitness over generations",
            FitnessMatplotlibChart(),
        )
        self.fitness_card.hide()
        self._content_layout.addWidget(self.fitness_card, 0, Qt.AlignTop | Qt.AlignLeft)

        self.diversity_card, self.diversity_chart, self.diversity_summary_label = self._make_card(
            "Population diversity (final generation)",
            TanimotoHistogramChart(),
        )
        self.diversity_card.hide()
        self._content_layout.addWidget(self.diversity_card, 0, Qt.AlignTop | Qt.AlignLeft)

        self._content_layout.addStretch(1)
        scroll.setWidget(content)

    def _make_card(self, title, chart_widget):
        card = QGroupBox(title)
        card.setObjectName("statsCard")
        card.setStyleSheet(STATS_CARD_QSS)
        card.setFixedWidth(self.CARD_WIDTH)
        card.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        _apply_card_shadow(card)

        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(8, 18, 8, 10)
        card_layout.setSpacing(8)
        card_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        card_layout.addWidget(chart_widget, 0, Qt.AlignLeft)

        summary_label = QLabel("")
        summary_label.setWordWrap(True)
        summary_label.setStyleSheet(
            "color: #424242; font-size: 12px; font-weight: normal; background: transparent;"
        )
        card_layout.addWidget(summary_label, 0, Qt.AlignLeft)
        return card, chart_widget, summary_label

    def update_from_statistics(self, stats: EvolutionStatistics):
        if not stats.has_data():
            self.empty_label.show()
            self.fitness_card.hide()
            self.diversity_card.hide()
            self.fitness_summary_label.setText("")
            self.diversity_summary_label.setText("")
            return

        self.empty_label.hide()
        self.fitness_card.show()
        self.fitness_chart.plot(
            stats.generation_indices,
            stats.best_fitness,
            stats.average_fitness,
        )
        self.fitness_summary_label.setText(
            fitness_summary_text(
                stats.generation_indices,
                stats.best_fitness,
                stats.average_fitness,
            )
        )

        final_entry = stats.final_generation_entry()
        diversity = final_entry.get("diversity") if final_entry else None
        if diversity and diversity.get("summary", {}).get("pair_count", 0) > 0:
            self.diversity_card.show()
            generation = final_entry.get("generation")
            self.diversity_chart.plot(
                diversity.get("pairwise_tanimoto", []),
                generation=generation,
            )
            self.diversity_summary_label.setText(
                tanimoto_summary_text(
                    diversity.get("summary"),
                    generation=generation,
                    similarities=diversity.get("pairwise_tanimoto", []),
                )
            )
        else:
            self.diversity_card.hide()
            self.diversity_summary_label.setText("")
