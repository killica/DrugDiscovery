from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor, QFont
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QGroupBox,
    QGraphicsDropShadowEffect,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QAbstractItemView,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from analysis.molecule_render import smiles_to_pixmap
from analysis.plots import (
    plot_fitness_over_generations,
    plot_tanimoto_histogram,
    fitness_summary_text,
    run_config_table_rows,
    operators_table_rows,
    best_molecule_summary_text,
    tanimoto_summary_text,
)
from evolutionRun import EvolutionStatistics
from moleculeBoxes import RoundedMoleculeImage


STATS_CARD_QSS = """
    QGroupBox#statsCard {
        background-color: rgb(255, 255, 255);
        border: 2px solid #43a047;
        border-radius: 15px;
        padding: 16px;
        padding-top: 28px;
        margin-top: 12px;
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

BODY_LABEL_QSS = (
    "color: #424242; font-size: 12px; font-weight: normal; "
    "background: transparent; padding: 2px 0;"
)

STATS_TABLE_QSS = """
    QTableWidget {
        border: 1px solid #e0e0e0;
        border-radius: 8px;
        background: #fafafa;
        font-size: 11px;
        color: #424242;
        gridline-color: #eeeeee;
    }
    QTableWidget::item {
        padding: 5px 8px;
    }
    QHeaderView::section {
        background: #e8f5e9;
        color: #1b5e20;
        font-weight: bold;
        border: none;
        border-bottom: 1px solid #c8e6c9;
        padding: 6px 8px;
    }
"""


def _apply_card_shadow(widget):
    shadow = QGraphicsDropShadowEffect(widget)
    shadow.setOffset(3, 5)
    shadow.setBlurRadius(20)
    shadow.setColor(QColor(0, 0, 0, 95))
    widget.setGraphicsEffect(shadow)


def _sync_wrapped_label(label):
    label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
    label.adjustSize()
    label.setMinimumHeight(label.sizeHint().height())


class _MatplotlibCardChart(FigureCanvas):
    CHART_WIDTH = 420
    CHART_HEIGHT = 240

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

    TOP_CARD_MIN_WIDTH = 260
    CHART_CARD_MIN_WIDTH = 400
    REPORT_MOLECULE_SIZE = 250

    def __init__(self):
        super().__init__()
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setStyleSheet("QScrollArea { background: transparent; border: none; }")
        outer.addWidget(scroll)

        content = QWidget()
        self._content_layout = QGridLayout(content)
        self._content_layout.setContentsMargins(16, 16, 16, 24)
        self._content_layout.setHorizontalSpacing(14)
        self._content_layout.setVerticalSpacing(16)
        self._content_layout.setColumnStretch(0, 1)
        self._content_layout.setColumnStretch(1, 1)
        self._content_layout.setColumnStretch(2, 1)

        self.empty_label = QLabel(
            "No evolution data yet.\nRun the genetic algorithm to collect statistics."
        )
        self.empty_label.setAlignment(Qt.AlignCenter)
        self.empty_label.setStyleSheet("color: #616161; font-size: 14px;")
        self._content_layout.addWidget(self.empty_label, 0, 0, 1, 3)

        self.config_card, self.config_table = self._make_table_card("Run settings", compact=True)
        self.config_card.hide()
        self._content_layout.addWidget(self.config_card, 0, 0, 1, 1, Qt.AlignTop)

        self.best_card, self.best_image, self.best_details_label = self._make_best_molecule_card()
        self.best_card.hide()
        self._content_layout.addWidget(self.best_card, 0, 1, 1, 1, Qt.AlignTop)

        self.operators_card, self.operators_table = self._make_table_card(
            "Genetic operators (full run)",
            compact=True,
        )
        self.operators_card.hide()
        self._content_layout.addWidget(self.operators_card, 0, 2, 1, 1, Qt.AlignTop)

        self.charts_row = QWidget()
        self.charts_layout = QHBoxLayout(self.charts_row)
        self.charts_layout.setContentsMargins(0, 0, 0, 0)
        self.charts_layout.setSpacing(14)

        self.fitness_card, self.fitness_chart, self.fitness_summary_label = self._make_card(
            "Fitness over generations",
            FitnessMatplotlibChart(),
        )
        self.fitness_card.hide()
        self.charts_layout.addWidget(self.fitness_card, 1)

        self.diversity_card, self.diversity_chart, self.diversity_summary_label = self._make_card(
            "Population diversity (final generation)",
            TanimotoHistogramChart(),
        )
        self.diversity_card.hide()
        self.charts_layout.addWidget(self.diversity_card, 1)

        self.charts_row.hide()
        self._content_layout.addWidget(self.charts_row, 1, 0, 1, 3)

        scroll.setWidget(content)

    def _apply_card_sizing(self, card, compact=False):
        min_width = self.TOP_CARD_MIN_WIDTH if compact else self.CHART_CARD_MIN_WIDTH
        card.setMinimumWidth(min_width)
        card.setMaximumWidth(16777215)
        card.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)

    def _make_body_label(self):
        label = QLabel("")
        label.setWordWrap(True)
        label.setStyleSheet(BODY_LABEL_QSS)
        return label

    def _make_card(self, title, chart_widget):
        card = QGroupBox(title)
        card.setObjectName("statsCard")
        card.setStyleSheet(STATS_CARD_QSS)
        self._apply_card_sizing(card)
        _apply_card_shadow(card)

        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(14, 12, 14, 14)
        card_layout.setSpacing(10)
        card_layout.setAlignment(Qt.AlignTop)
        card_layout.addWidget(chart_widget, 0, Qt.AlignHCenter)

        summary_label = self._make_body_label()
        card_layout.addWidget(summary_label, 0, Qt.AlignLeft)
        return card, chart_widget, summary_label

    def _make_stats_table(self):
        table = QTableWidget(0, 2)
        table.setHorizontalHeaderLabels(["Parameter", "Value"])
        table.horizontalHeader().setStretchLastSection(True)
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        table.verticalHeader().setVisible(False)
        table.setShowGrid(True)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionMode(QAbstractItemView.NoSelection)
        table.setFocusPolicy(Qt.NoFocus)
        table.setStyleSheet(STATS_TABLE_QSS)
        table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        return table

    def _make_table_card(self, title, compact=False):
        card = QGroupBox(title)
        card.setObjectName("statsCard")
        card.setStyleSheet(STATS_CARD_QSS)
        self._apply_card_sizing(card, compact=compact)
        _apply_card_shadow(card)

        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(14, 12, 14, 14)
        card_layout.setSpacing(6)
        card_layout.setAlignment(Qt.AlignTop)

        table = self._make_stats_table()
        card_layout.addWidget(table, 0, Qt.AlignLeft)
        return card, table

    def _populate_stats_table(self, table, rows, empty_message=None):
        table.clearContents()
        table.setRowCount(0)

        if not rows:
            table.setColumnCount(2)
            table.setHorizontalHeaderLabels(["Parameter", "Value"])
            if empty_message:
                table.setRowCount(1)
                message_item = QTableWidgetItem(empty_message)
                message_item.setFlags(Qt.ItemIsEnabled)
                table.setItem(0, 0, message_item)
                table.setSpan(0, 0, 1, 2)
                table.resizeRowsToContents()
                table.setFixedHeight(table.horizontalHeader().height() + table.rowHeight(0) + 6)
            else:
                table.setRowCount(0)
                table.setFixedHeight(table.horizontalHeader().height() + 6)
            return

        table.setColumnCount(2)
        table.setHorizontalHeaderLabels(["Parameter", "Value"])
        table.horizontalHeader().setStretchLastSection(True)
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)

        section_font = QFont()
        section_font.setBold(True)

        for label, value in rows:
            row = table.rowCount()
            table.insertRow(row)

            if label == "__section__":
                section_item = QTableWidgetItem(str(value))
                section_item.setFlags(Qt.ItemIsEnabled)
                section_item.setFont(section_font)
                section_item.setBackground(QColor("#e8f5e9"))
                section_item.setForeground(QColor("#1b5e20"))
                table.setItem(row, 0, section_item)
                table.setSpan(row, 0, 1, 2)
                continue

            key_item = QTableWidgetItem(str(label))
            key_item.setFlags(Qt.ItemIsEnabled)
            value_item = QTableWidgetItem(str(value))
            value_item.setFlags(Qt.ItemIsEnabled)
            value_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            table.setItem(row, 0, key_item)
            table.setItem(row, 1, value_item)

        table.resizeRowsToContents()
        table.resizeColumnToContents(0)
        total_height = table.horizontalHeader().height() + sum(
            table.rowHeight(row_index) for row_index in range(table.rowCount())
        ) + 6
        table.setFixedHeight(total_height)

    def _make_best_molecule_card(self):
        card = QGroupBox("Best molecule so far")
        card.setObjectName("statsCard")
        card.setStyleSheet(STATS_CARD_QSS)
        self._apply_card_sizing(card, compact=True)
        _apply_card_shadow(card)

        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(10, 10, 10, 10)
        card_layout.setSpacing(8)
        card_layout.setAlignment(Qt.AlignTop | Qt.AlignHCenter)

        image_frame = QWidget()
        image_frame.setFixedSize(self.REPORT_MOLECULE_SIZE, self.REPORT_MOLECULE_SIZE)
        image_layout = QVBoxLayout(image_frame)
        image_layout.setContentsMargins(0, 0, 0, 0)

        image_label = RoundedMoleculeImage(None)
        image_label.setFixedSize(self.REPORT_MOLECULE_SIZE, self.REPORT_MOLECULE_SIZE)
        image_layout.addWidget(image_label)

        image_shadow = QGraphicsDropShadowEffect(image_frame)
        image_shadow.setOffset(3, 5)
        image_shadow.setBlurRadius(20)
        image_shadow.setColor(QColor(0, 0, 0, 95))
        image_frame.setGraphicsEffect(image_shadow)

        image_row = QHBoxLayout()
        image_row.addStretch(1)
        image_row.addWidget(image_frame, 0, Qt.AlignHCenter)
        image_row.addStretch(1)
        card_layout.addLayout(image_row)

        details_label = self._make_body_label()
        details_label.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
        details_label.setStyleSheet(
            "color: #424242; font-size: 11px; font-weight: normal; "
            "background: transparent; padding: 2px 0;"
        )
        card_layout.addWidget(details_label, 0, Qt.AlignHCenter)

        return card, image_label, details_label

    def _set_card_text(self, label, text):
        label.setText(text)
        _sync_wrapped_label(label)

    def _resolve_best_molecule(self, stats, best_molecule):
        if best_molecule:
            return best_molecule
        return getattr(stats, "best_candidate", None)

    def update_from_statistics(self, stats: EvolutionStatistics, best_molecule=None):
        if not stats.has_data():
            self.empty_label.show()
            self.best_card.hide()
            self.config_card.hide()
            self.charts_row.hide()
            self.fitness_card.hide()
            self.operators_card.hide()
            self.diversity_card.hide()
            self._populate_stats_table(self.config_table, [])
            self.fitness_summary_label.setText("")
            self._populate_stats_table(self.operators_table, [])
            self.diversity_summary_label.setText("")
            self.best_details_label.setText("")
            return

        self.empty_label.hide()

        candidate = self._resolve_best_molecule(stats, best_molecule)
        if candidate and candidate.get("smiles"):
            self.best_card.show()
            pixmap = smiles_to_pixmap(candidate["smiles"], size=self.REPORT_MOLECULE_SIZE)
            self.best_image.setMoleculePixmap(pixmap)
            self._set_card_text(
                self.best_details_label,
                best_molecule_summary_text(candidate),
            )
        else:
            self.best_card.hide()
            self.best_details_label.setText("")

        initial_population_size = None
        if stats.generations_data:
            initial_population_size = stats.generations_data[0].get("population_size")

        self.config_card.show()
        self._populate_stats_table(
            self.config_table,
            run_config_table_rows(stats.config, initial_population_size=initial_population_size),
        )

        self.charts_row.show()
        self.fitness_card.show()
        self.fitness_chart.plot(
            stats.generation_indices,
            stats.best_fitness,
            stats.average_fitness,
        )
        self._set_card_text(
            self.fitness_summary_label,
            fitness_summary_text(
                stats.generation_indices,
                stats.best_fitness,
                stats.average_fitness,
            ),
        )

        self.operators_card.show()
        self._populate_stats_table(
            self.operators_table,
            operators_table_rows(stats.operators_summary),
            empty_message=(
                "Operator statistics not available yet. "
                "Complete at least one full GA run (Launch or Generate)."
            ),
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
            self._set_card_text(
                self.diversity_summary_label,
                tanimoto_summary_text(
                    diversity.get("summary"),
                    generation=generation,
                    similarities=diversity.get("pairwise_tanimoto", []),
                ),
            )
        else:
            self.diversity_card.hide()
            self.diversity_summary_label.setText("")
