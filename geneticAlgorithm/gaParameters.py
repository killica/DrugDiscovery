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
)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon
import geneticAlgorithm
from GAConfig import CrossoverMode, MutationMode
from mutationInfo import MutationInfo

EVOLUTION_ACTION_BTN_WIDTH = 200

GA_PROGRESS_BAR_STYLE = """
QProgressBar {
    border: 1px solid #9e9e9e;
    border-radius: 6px;
    background-color: #e8e8e8;
    text-align: center;
    min-height: 24px;
    max-height: 28px;
    color: #1a1a1a;
}
QProgressBar::chunk {
    background-color: #2e7d32;
    border-radius: 5px;
}
"""


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

    def showConfiguredParametersDialog(self, parent=None):
        cfg = self.application.gaConfig
        if cfg.rouletteSelection:
            selection = "Roulette wheel"
        else:
            selection = f"Tournament (size {cfg.tournamentSize})"
        crossover = cfg.crossoverMode.name
        mutation = cfg.mutationMode.name
        QMessageBox.information(
            parent or self.application,
            "Genetic algorithm parameters",
            (
                f"Number of generations: {cfg.generations}\n"
                f"Selection: {selection}\n"
                f"Elitism size: {cfg.elitismSize}\n"
                f"Mutation probability: {cfg.mutationProbability}\n"
                f"Crossover representation: {crossover}\n"
                f"Mutation representation: {mutation}"
            ),
        )

    def onLaunchButtonClicked(self):
        # moleculeBoxes is a reference to the right half of the scene
        moleculeBoxes = self.application.moleculeBoxes

        moleculeBoxes.generateButton = QPushButton("Generate next", moleculeBoxes.application)
        moleculeBoxes.generateButton.setFixedWidth(EVOLUTION_ACTION_BTN_WIDTH)
        moleculeBoxes.generateButton.clicked.connect(moleculeBoxes.onGenerateButtonClicked)
        moleculeBoxes.generateButton.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)

        moleculeBoxes.finalButton = QPushButton("Jump to final")
        moleculeBoxes.finalButton.setFixedWidth(EVOLUTION_ACTION_BTN_WIDTH)
        moleculeBoxes.finalButton.clicked.connect(moleculeBoxes.onFinalButtonClicked)
        moleculeBoxes.finalButton.setStyleSheet("""
            QPushButton {
                background-color: #6495ED;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #0047AB;
            }
        """)

        moleculeBoxes.saveButton = QPushButton("Save the best")
        moleculeBoxes.saveButton.setFixedWidth(EVOLUTION_ACTION_BTN_WIDTH)
        moleculeBoxes.saveButton.clicked.connect(moleculeBoxes.onSaveButtonClicked)
        moleculeBoxes.saveButton.setStyleSheet("""
            QPushButton {
                background-color: #606060;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: black;
            }
        """)

        moleculeBoxes.saveLabel = QLabel("Saved!")
        moleculeBoxes.saveLabel.setStyleSheet("color: transparent; font-style: italic;")

        moleculeBoxes.saveBox = QHBoxLayout()
        moleculeBoxes.saveBox.setContentsMargins(0, 0, 0, 0)
        moleculeBoxes.saveBox.setSpacing(8)
        moleculeBoxes.saveBox.addWidget(moleculeBoxes.saveButton, 0, Qt.AlignLeft)
        moleculeBoxes.saveBox.addWidget(moleculeBoxes.saveLabel, 0, Qt.AlignVCenter)
        moleculeBoxes.saveBox.addStretch(1)

        moleculeBoxes.saveCnt = QWidget()
        moleculeBoxes.saveCnt.setLayout(moleculeBoxes.saveBox)
        moleculeBoxes.saveCnt.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        moleculeBoxes.restartButton = QPushButton("Restart analysis")
        moleculeBoxes.restartButton.setFixedWidth(EVOLUTION_ACTION_BTN_WIDTH)
        moleculeBoxes.restartButton.clicked.connect(moleculeBoxes.onRestartButtonClicked)
        moleculeBoxes.restartButton.setStyleSheet("""
            QPushButton {
                background-color: #ff4040;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: red;
            }
        """)

        moleculeBoxes.gaParametersButton = QPushButton("GA parameters")
        moleculeBoxes.gaParametersButton.setFixedWidth(EVOLUTION_ACTION_BTN_WIDTH)
        moleculeBoxes.gaParametersButton.clicked.connect(
            lambda: self.showConfiguredParametersDialog(moleculeBoxes.application)
        )
        moleculeBoxes.gaParametersButton.setStyleSheet("""
            QPushButton {
                background-color: #455a64;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #37474f;
            }
        """)

        moleculeBoxes.rightVBox3 = QVBoxLayout()
        moleculeBoxes.rightVBox3.setSpacing(10)
        moleculeBoxes.rightVBox3.setContentsMargins(0, 0, 0, 0)
        moleculeBoxes.rightVBox3.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        for w in (
            moleculeBoxes.generateButton,
            moleculeBoxes.finalButton,
            moleculeBoxes.saveCnt,
            moleculeBoxes.restartButton,
        ):
            moleculeBoxes.rightVBox3.addWidget(w, 0, Qt.AlignLeft)

        moleculeBoxes.rightBtnCnt = QWidget()
        moleculeBoxes.rightBtnCnt.setLayout(moleculeBoxes.rightVBox3)
        moleculeBoxes.rightBtnCnt.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Maximum)

        moleculeBoxes.bestBox = moleculeBoxes.createMoleculeBox("", "To be determined", 0.0, 0, -1)

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

        moleculeBoxes.progressVBox = QVBoxLayout()

        moleculeBoxes.generationLabel = QLabel(f"Generation: 2/{self.application.gaConfig.generations}")
        moleculeBoxes.generationLabel.setStyleSheet("color: blue; font-style: bold;")

        moleculeBoxes.generationProgress = QProgressBar()
        moleculeBoxes.generationProgress.setRange(0, self.application.gaConfig.generations)
        moleculeBoxes.generationProgress.setValue(2)
        moleculeBoxes.generationProgress.setStyleSheet(GA_PROGRESS_BAR_STYLE)
        moleculeBoxes.generationProgress.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        moleculeBoxes.generationProgress.setTextVisible(True)

        n_selected = len(moleculeBoxes.selectedMolecules)
        moleculeBoxes.individualLabel = QLabel(f"Individual: 0/{n_selected}")
        moleculeBoxes.individualLabel.setStyleSheet("color: blue; font-style: bold;")

        moleculeBoxes.individualProgress = QProgressBar()
        moleculeBoxes.individualProgress.setRange(0, max(n_selected, 1))
        moleculeBoxes.individualProgress.setValue(0)
        moleculeBoxes.individualProgress.setStyleSheet(GA_PROGRESS_BAR_STYLE)
        moleculeBoxes.individualProgress.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        moleculeBoxes.individualProgress.setTextVisible(True)

        moleculeBoxes.progressVBox.addWidget(moleculeBoxes.generationLabel)
        moleculeBoxes.progressVBox.addWidget(moleculeBoxes.generationProgress)
        moleculeBoxes.progressVBox.addSpacing(30)
        moleculeBoxes.progressVBox.addWidget(moleculeBoxes.individualLabel)
        moleculeBoxes.progressVBox.addWidget(moleculeBoxes.individualProgress)

        moleculeBoxes.progressCnt = QWidget()
        moleculeBoxes.progressCnt.setLayout(moleculeBoxes.progressVBox)
        moleculeBoxes.progressCnt.setMinimumWidth(300)
        moleculeBoxes.progressCnt.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

        moleculeBoxes.evolutionControlsLayout.addWidget(
            moleculeBoxes.rightBtnCnt, 0, Qt.AlignTop | Qt.AlignLeft
        )
        moleculeBoxes.evolutionControlsLayout.addWidget(
            moleculeBoxes.bestBox, 0, Qt.AlignTop | Qt.AlignLeft
        )
        moleculeBoxes.evolutionControlsLayout.addWidget(
            moleculeBoxes.progressCnt, 0, Qt.AlignTop | Qt.AlignLeft
        )
        moleculeBoxes.evolutionControlsLayout.addStretch(1)
        moleculeBoxes.evolutionControlsLayout.addWidget(
            moleculeBoxes.gaParametersButton, 0, Qt.AlignLeft
        )

        moleculeBoxes.rightCont3.setLayout(moleculeBoxes.evolutionControlsLayout)
        moleculeBoxes.rightCont3.setMinimumWidth(380)
        moleculeBoxes.rightCont3.setMaximumWidth(520)
        moleculeBoxes.rightCont3.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

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
            )
        finally:
            moleculeBoxes._ga_running = False
            moleculeBoxes._set_evolution_actions_enabled(True)

        if getattr(self.application, "_cancel_evolution", False):
            return

        moleculeBoxes.loadNewGeneration(tuple(self.application.sliderValues))
        moleculeBoxes._mark_evolution_complete_if_done()
        
        

