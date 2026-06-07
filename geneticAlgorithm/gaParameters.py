from PyQt5.QtWidgets import (
    QApplication,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QLabel,
    QSpinBox,
    QWidget,
    QPushButton,
    QLineEdit,
    QCheckBox,
    QRadioButton,
    QButtonGroup,
    QProgressBar,
    QSizePolicy,
    QMessageBox,
    QGroupBox,
    QGraphicsDropShadowEffect,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
import geneticAlgorithm
from GAConfig import CrossoverMode, MutationMode
from mutationInfo import MutationInfo

EVOLUTION_PANEL_QSS = """
    QGroupBox#evolutionPanel {
        background-color: rgb(255, 255, 255);
        border: 2px solid #43a047;
        border-radius: 12px;
        padding: 12px;
        padding-top: 24px;
        margin-top: 8px;
        font-size: 13px;
        font-weight: bold;
        color: #1b5e20;
    }
    QGroupBox#evolutionPanel::title {
        subcontrol-origin: margin;
        left: 12px;
        padding: 0 6px;
    }
"""

EVOLUTION_PANEL_WIDTH = 256


def _configure_evolution_panel_card(card):
    card.setObjectName("evolutionPanel")
    card.setStyleSheet(EVOLUTION_PANEL_QSS)
    card.setFixedWidth(EVOLUTION_PANEL_WIDTH)
    card.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
    _apply_evolution_panel_shadow(card)


GA_PROGRESS_BAR_STYLE = """
QProgressBar {
    border: 1px solid #c5cae9;
    border-radius: 6px;
    background-color: #f5f5f5;
    text-align: center;
    min-height: 22px;
    max-height: 26px;
    color: #1a1a1a;
    font-size: 11px;
}
QProgressBar::chunk {
    background-color: #2e7d32;
    border-radius: 5px;
}
"""

GA_PROGRESS_INDIVIDUAL_STYLE = """
QProgressBar {
    border: 1px solid #c5cae9;
    border-radius: 6px;
    background-color: #f5f5f5;
    text-align: center;
    min-height: 22px;
    max-height: 26px;
    color: #1a1a1a;
    font-size: 11px;
}
QProgressBar::chunk {
    background-color: #1565c0;
    border-radius: 5px;
}
"""


def _evolution_action_button_style(background, hover):
    return f"""
        QPushButton {{
            background-color: {background};
            color: white;
            border: none;
            border-radius: 8px;
            padding: 5px 10px;
            font-size: 12px;
            font-weight: bold;
            min-height: 28px;
            max-height: 32px;
        }}
        QPushButton:hover {{
            background-color: {hover};
        }}
        QPushButton:disabled {{
            background-color: #bdbdbd;
            color: #f5f5f5;
        }}
    """


def _apply_evolution_panel_shadow(widget):
    shadow = QGraphicsDropShadowEffect(widget)
    shadow.setOffset(2, 4)
    shadow.setBlurRadius(16)
    shadow.setColor(QColor(0, 0, 0, 70))
    widget.setGraphicsEffect(shadow)


class GAParameters:
    def __init__(self, application):
        self.application = application
        self.paramsLayout = QVBoxLayout()
        self.paramsLayout.setSpacing(12)
        self.paramsLayout.setContentsMargins(0, 0, 0, 0)

        self.numGenLabel = QLabel("Number of generations:", application)
        self.tournamentSizeLabel = QLabel("Tournament size:", application)
        self.elitismSizeLabel = QLabel("Elitism size:", application)
        self.mutationProbabilityLabel = QLabel("Mutation probability:", application)
        self.crossoverLabel = QLabel("Crossover:", application)
        self.mutationLabel = QLabel("Mutation:", application)
        self.generationSpin = QSpinBox(application)
        self.generationSpin.setMaximum(200)
        self.generationSpin.setFixedWidth(70)
        self.generationSpin.setValue(100)
        self.rouletteCheckBox = QCheckBox("Roulette selection", application)
        self.rouletteCheckBox.setStyleSheet("""
            QCheckBox {
                text-decoration: none; 
            }
            QCheckBox::indicator {
                width: 20px;
                height: 20px;
                border: 2px solid #777;
                border-radius: 5px;
                background-color: white;
            }
            QCheckBox::indicator:checked {
                background-color: lightgray; 
                border: 2px solid green;
            }
            QCheckBox::indicator:unchecked {
                background-color: white;
                border: 2px solid #777;
            }
        """)

        self.rouletteCheckBox.stateChanged.connect(self.onRouletteChanged)

        for label in (
            self.numGenLabel,
            self.tournamentSizeLabel,
            self.elitismSizeLabel,
            self.mutationProbabilityLabel,
            self.crossoverLabel,
            self.mutationLabel,
        ):
            label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

        self.tournamentSpin = QSpinBox(application)
        self.tournamentSpin.setMaximum(100)
        self.tournamentSpin.setFixedWidth(70)
        self.tournamentSpin.setValue(7)
        self.elitismSpin = QSpinBox(application)
        self.elitismSpin.setMaximum(100)
        self.elitismSpin.setFixedWidth(70)
        self.elitismSpin.setValue(2)
        self.mutationLineEdit = QLineEdit(application)
        self.mutationLineEdit.setFixedWidth(70)
        self.mutationLineEdit.setText("0.05")

        tournament_field = QWidget(application)
        tournament_field.setMinimumWidth(280)
        tournament_row = QHBoxLayout(tournament_field)
        tournament_row.setContentsMargins(0, 0, 0, 0)
        tournament_row.setSpacing(16)
        tournament_row.addWidget(self.tournamentSpin, 0, Qt.AlignLeft)
        tournament_row.addWidget(self.rouletteCheckBox, 0, Qt.AlignLeft)
        tournament_row.addStretch(1)

        self.formLayout = QFormLayout()
        self.formLayout.setSpacing(10)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.formLayout.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.formLayout.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.formLayout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        self.formLayout.setRowWrapPolicy(QFormLayout.DontWrapRows)
        self.formLayout.addRow(self.numGenLabel, self.generationSpin)
        self.formLayout.addRow(self.tournamentSizeLabel, tournament_field)
        self.formLayout.addRow(self.elitismSizeLabel, self.elitismSpin)
        self.formLayout.addRow(self.mutationProbabilityLabel, self.mutationLineEdit)

        self.crossoverSmilesRadio = QRadioButton("SMILES", application)
        self.crossoverSelfiesRadio = QRadioButton("SELFIES", application)
        self.crossoverGraphRadio = QRadioButton("Graph", application)
        self.crossoverBricsRadio = QRadioButton("BRICS", application)
        self.crossoverSelfiesRadio.setChecked(True)
        self.crossoverGroup = QButtonGroup(application)
        self.crossoverGroup.addButton(self.crossoverSmilesRadio, CrossoverMode.SMILES.value)
        self.crossoverGroup.addButton(self.crossoverSelfiesRadio, CrossoverMode.SELFIES.value)
        self.crossoverGroup.addButton(self.crossoverGraphRadio, CrossoverMode.GRAPH.value)
        self.crossoverGroup.addButton(self.crossoverBricsRadio, CrossoverMode.BRICS.value)

        crossover_field = QWidget(application)
        crossover_field.setMinimumWidth(360)
        crossover_row = QHBoxLayout(crossover_field)
        crossover_row.setContentsMargins(0, 0, 0, 0)
        crossover_row.setSpacing(16)
        crossover_row.addWidget(self.crossoverSmilesRadio, 0, Qt.AlignLeft)
        crossover_row.addWidget(self.crossoverSelfiesRadio, 0, Qt.AlignLeft)
        crossover_row.addWidget(self.crossoverGraphRadio, 0, Qt.AlignLeft)
        crossover_row.addWidget(self.crossoverBricsRadio, 0, Qt.AlignLeft)
        crossover_row.addStretch(1)
        self.formLayout.addRow(self.crossoverLabel, crossover_field)

        self.mutationSmilesRadio = QRadioButton("SMILES", application)
        self.mutationSelfiesRadio = QRadioButton("SELFIES", application)
        self.mutationGraphRadio = QRadioButton("Graph", application)
        self.mutationBricsRadio = QRadioButton("BRICS", application)
        self.mutationSmilesRadio.setChecked(True)
        self.mutationGroup = QButtonGroup(application)
        self.mutationGroup.addButton(self.mutationSmilesRadio, MutationMode.SMILES.value)
        self.mutationGroup.addButton(self.mutationSelfiesRadio, MutationMode.SELFIES.value)
        self.mutationGroup.addButton(self.mutationGraphRadio, MutationMode.GRAPH.value)
        self.mutationGroup.addButton(self.mutationBricsRadio, MutationMode.BRICS.value)

        mutation_field = QWidget(application)
        mutation_field.setMinimumWidth(360)
        mutation_row = QHBoxLayout(mutation_field)
        mutation_row.setContentsMargins(0, 0, 0, 0)
        mutation_row.setSpacing(16)
        mutation_row.addWidget(self.mutationSmilesRadio, 0, Qt.AlignLeft)
        mutation_row.addWidget(self.mutationSelfiesRadio, 0, Qt.AlignLeft)
        mutation_row.addWidget(self.mutationGraphRadio, 0, Qt.AlignLeft)
        mutation_row.addWidget(self.mutationBricsRadio, 0, Qt.AlignLeft)
        mutation_row.addStretch(1)
        self.formLayout.addRow(self.mutationLabel, mutation_field)

        self.paramsLayout.addLayout(self.formLayout)

        self.launchButton = QPushButton("Launch search!", application)
        self.launchButton.setFixedWidth(200)
        self.launchButton.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 16px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)
        self.launchButton.clicked.connect(self.onLaunchButtonClicked)

        launch_row = QHBoxLayout()
        launch_row.setContentsMargins(0, 8, 0, 0)
        launch_row.addWidget(self.launchButton, 0, Qt.AlignLeft)
        launch_row.addStretch(1)
        self.paramsLayout.addLayout(launch_row)

        self.container = QWidget()
        self.container.setLayout(self.paramsLayout)
        self.container.setMinimumWidth(460)
        self.container.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)

    def onRouletteChanged(self, state):
        if state == 2:
            self.tournamentSpin.setDisabled(True)
        else:
            self.tournamentSpin.setDisabled(False)

    def getGAParametersWidget(self):
        return self.container

    def selectedCrossoverMode(self) -> CrossoverMode:
        return CrossoverMode(self.crossoverGroup.checkedId())

    def selectedMutationMode(self) -> MutationMode:
        return MutationMode(self.mutationGroup.checkedId())

    def onLaunchButtonClicked(self):
        # moleculeBoxes is a reference to the right half of the scene
        moleculeBoxes = self.application.moleculeBoxes

        moleculeBoxes.generateButton = QPushButton("Generate next", moleculeBoxes.application)
        moleculeBoxes.generateButton.clicked.connect(moleculeBoxes.onGenerateButtonClicked)
        moleculeBoxes.generateButton.setStyleSheet(
            _evolution_action_button_style("#43a047", "#2e7d32")
        )

        moleculeBoxes.finalButton = QPushButton("Jump to final")
        moleculeBoxes.finalButton.clicked.connect(moleculeBoxes.onFinalButtonClicked)
        moleculeBoxes.finalButton.setStyleSheet(
            _evolution_action_button_style("#5c6bc0", "#3949ab")
        )

        moleculeBoxes.showStatsButton = QPushButton("Show report")
        moleculeBoxes.showStatsButton.clicked.connect(moleculeBoxes.onShowStatsButtonClicked)
        moleculeBoxes.showStatsButton.setStyleSheet(
            _evolution_action_button_style("#455a64", "#37474f")
        )

        moleculeBoxes.restartButton = QPushButton("Restart analysis")
        moleculeBoxes.restartButton.clicked.connect(moleculeBoxes.onRestartButtonClicked)
        moleculeBoxes.restartButton.setStyleSheet(
            _evolution_action_button_style("#e53935", "#c62828")
        )

        action_buttons = (
            moleculeBoxes.generateButton,
            moleculeBoxes.finalButton,
            moleculeBoxes.showStatsButton,
            moleculeBoxes.restartButton,
        )
        for button in action_buttons:
            button.setMinimumWidth(220)
            button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        actions_card = QGroupBox("Actions")
        _configure_evolution_panel_card(actions_card)

        actions_layout = QVBoxLayout(actions_card)
        actions_layout.setContentsMargins(12, 10, 12, 12)
        actions_layout.setSpacing(8)
        for button in action_buttons:
            actions_layout.addWidget(button, 0, Qt.AlignTop)

        moleculeBoxes.actionsCard = actions_card

        self.launchButton.setDisabled(True)
        self.launchButton.setStyleSheet("""
            QPushButton {
                background-color: #757575;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 16px;
                font-weight: bold;
            }
        """)

        # Necessary parameters for genetic algorithm
        self.application.gaConfig.rouletteSelection = self.rouletteCheckBox.isChecked()
        self.application.gaConfig.generations = self.generationSpin.value()
        self.application.gaConfig.tournamentSize = self.tournamentSpin.value()
        self.application.gaConfig.elitismSize = self.elitismSpin.value()
        self.application.gaConfig.mutationProbability = float(self.mutationLineEdit.text())
        self.application.gaConfig.crossoverMode = self.selectedCrossoverMode()
        self.application.gaConfig.mutationMode = self.selectedMutationMode()

        progress_card = QGroupBox("Progress")
        _configure_evolution_panel_card(progress_card)

        progress_layout = QVBoxLayout(progress_card)
        progress_layout.setContentsMargins(12, 10, 12, 12)
        progress_layout.setSpacing(8)

        moleculeBoxes.generationLabel = QLabel(
            f"Generation: 2/{self.application.gaConfig.generations}"
        )
        moleculeBoxes.generationLabel.setStyleSheet(
            "color: #1b5e20; font-size: 12px; font-weight: bold;"
        )

        moleculeBoxes.generationProgress = QProgressBar()
        moleculeBoxes.generationProgress.setRange(0, self.application.gaConfig.generations)
        moleculeBoxes.generationProgress.setValue(2)
        moleculeBoxes.generationProgress.setStyleSheet(GA_PROGRESS_BAR_STYLE)
        moleculeBoxes.generationProgress.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        moleculeBoxes.generationProgress.setTextVisible(True)

        n_selected = len(moleculeBoxes.selectedMolecules)
        moleculeBoxes.individualLabel = QLabel(f"Individual: 0/{n_selected}")
        moleculeBoxes.individualLabel.setStyleSheet(
            "color: #1565c0; font-size: 12px; font-weight: bold;"
        )

        moleculeBoxes.individualProgress = QProgressBar()
        moleculeBoxes.individualProgress.setRange(0, max(n_selected, 1))
        moleculeBoxes.individualProgress.setValue(0)
        moleculeBoxes.individualProgress.setStyleSheet(GA_PROGRESS_INDIVIDUAL_STYLE)
        moleculeBoxes.individualProgress.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        moleculeBoxes.individualProgress.setTextVisible(True)

        progress_layout.addWidget(moleculeBoxes.generationLabel)
        progress_layout.addWidget(moleculeBoxes.generationProgress)
        progress_layout.addSpacing(6)
        progress_layout.addWidget(moleculeBoxes.individualLabel)
        progress_layout.addWidget(moleculeBoxes.individualProgress)

        moleculeBoxes.progressCard = progress_card

        moleculeBoxes.evolutionControlsLayout.setSpacing(14)
        moleculeBoxes.evolutionControlsLayout.setAlignment(Qt.AlignCenter)
        moleculeBoxes.evolutionControlsLayout.addStretch(1)
        moleculeBoxes.evolutionControlsLayout.addWidget(
            moleculeBoxes.actionsCard, 0, Qt.AlignCenter
        )
        moleculeBoxes.evolutionControlsLayout.addWidget(
            moleculeBoxes.progressCard, 0, Qt.AlignCenter
        )
        moleculeBoxes.evolutionControlsLayout.addStretch(1)

        moleculeBoxes.rightCont3.setLayout(moleculeBoxes.evolutionControlsLayout)
        moleculeBoxes.rightCont3.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Show evolution / progress on stage 3 before the GA blocks the event loop.
        self.application.show_stage_3()

        self.application.sbmtBtn.setDisabled(True)
        self.application.resBtn.setDisabled(True)

        self.application.blockTransfer = True
        self.application._cancel_evolution = False

        # Paint progress widgets once before the long-running GA blocks the event loop.
        QApplication.processEvents()

        moleculeBoxes._ga_running = True
        moleculeBoxes._set_evolution_actions_enabled(False)
        geneticAlgorithm.reset_crossover_stats()
        geneticAlgorithm.reset_mutation_stats()
        cfg = self.application.gaConfig
        self.application.evolution_statistics.begin_run(
            {
                "generations": cfg.generations,
                "tournament_size": cfg.tournamentSize,
                "elitism_size": cfg.elitismSize,
                "mutation_probability": cfg.mutationProbability,
                "roulette_selection": cfg.rouletteSelection,
                "crossover_mode": cfg.crossoverMode.name,
                "mutation_mode": cfg.mutationMode.name,
            }
        )
        try:
            moleculeBoxes.newGenerationMolecules = geneticAlgorithm.geneticAlgorithm(
                moleculeBoxes.selectedMolecules,
                True,
                self.application.gaConfig.generations,
                self.application.gaConfig.rouletteSelection,
                self.application.gaConfig.tournamentSize,
                self.application.gaConfig.elitismSize,
                self.application.gaConfig.mutationProbability,
                self.application.gaConfig.crossoverMode,
                self.application.gaConfig.mutationMode,
                self.application.getMutationInfo(),
                moleculeBoxes.individualLabel,
                moleculeBoxes.individualProgress,
                cancel_check=lambda: getattr(self.application, "_cancel_evolution", False),
                on_generation_start=moleculeBoxes._live_generation_on_start,
                on_new_individual=moleculeBoxes._live_generation_on_new_individual,
                evolution_stats=self.application.evolution_statistics,
                record_initial=True,
            )
        finally:
            moleculeBoxes._ga_running = False
            moleculeBoxes._set_evolution_actions_enabled(True)

        if getattr(self.application, "_cancel_evolution", False):
            return

        moleculeBoxes.loadNewGeneration(tuple(self.application.sliderValues))
        moleculeBoxes._mark_evolution_complete_if_done()
        
        

