from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class EvolutionStatistics:
    """Stores metrics collected during a genetic algorithm run."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.generation_indices = []
        self.best_qed = []
        self.average_qed = []
        self._next_generation = 1

    def has_data(self):
        return len(self.generation_indices) > 0

    def record_generation(self, population):
        if not population:
            return
        qeds = [individual.getQED() for individual in population]
        self.generation_indices.append(self._next_generation)
        self.best_qed.append(max(qeds))
        self.average_qed.append(sum(qeds) / len(qeds))
        self._next_generation += 1


class QedMatplotlibChart(FigureCanvas):
    """Compact matplotlib chart for best and average QED over generations."""

    def __init__(self, parent=None):
        self._figure = Figure(figsize=(5.0, 2.2), dpi=100)
        self._axes = self._figure.add_subplot(111)
        super().__init__(self._figure)
        self.setParent(parent)
        self.setFixedHeight(240)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._figure.tight_layout(pad=1.2)

    def plot(self, generations, best_qed, average_qed):
        self._axes.clear()
        self._axes.plot(
            generations,
            best_qed,
            marker="o",
            markersize=4,
            linewidth=1.6,
            color="#2e7d32",
            label="Best QED",
        )
        self._axes.plot(
            generations,
            average_qed,
            marker="o",
            markersize=4,
            linewidth=1.6,
            color="#1565c0",
            label="Average QED",
        )
        self._axes.set_xlabel("Generation", fontsize=9)
        self._axes.set_ylabel("QED", fontsize=9)
        self._axes.tick_params(labelsize=8)
        self._axes.grid(True, linestyle="--", alpha=0.35)
        self._axes.legend(loc="best", fontsize=8)
        self._figure.tight_layout(pad=1.2)
        self.draw()


class EvolutionStatsChart(QWidget):
    """Scrollable report panel; add more chart sections here over time."""

    def __init__(self):
        super().__init__()
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(14)

        self.empty_label = QLabel(
            "No evolution data yet.\nRun the genetic algorithm to collect statistics."
        )
        self.empty_label.setAlignment(Qt.AlignCenter)
        self.empty_label.setStyleSheet("color: #616161; font-size: 14px;")
        layout.addWidget(self.empty_label)

        self.qed_section = QWidget()
        qed_layout = QVBoxLayout(self.qed_section)
        qed_layout.setContentsMargins(0, 0, 0, 0)
        qed_layout.setSpacing(6)

        qed_title = QLabel("QED over generations")
        qed_title.setStyleSheet("font-size: 14px; font-weight: bold; color: #37474f;")
        qed_layout.addWidget(qed_title)

        self.qed_chart = QedMatplotlibChart()
        qed_layout.addWidget(self.qed_chart)

        self.summary_label = QLabel("")
        self.summary_label.setStyleSheet("color: #424242; font-size: 12px;")
        qed_layout.addWidget(self.summary_label)

        self.qed_section.hide()
        layout.addWidget(self.qed_section, 0, Qt.AlignTop)
        layout.addStretch(1)

        scroll.setWidget(content)

    def update_from_statistics(self, stats: EvolutionStatistics):
        if not stats.has_data():
            self.empty_label.show()
            self.qed_section.hide()
            self.summary_label.setText("")
            return

        self.empty_label.hide()
        self.qed_section.show()
        self.qed_chart.plot(
            stats.generation_indices,
            stats.best_qed,
            stats.average_qed,
        )
        self.summary_label.setText(
            f"Generations recorded: {len(stats.generation_indices)}  |  "
            f"Best QED: {max(stats.best_qed):.4f}  |  "
            f"Final average QED: {stats.average_qed[-1]:.4f}"
        )
