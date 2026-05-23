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
)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon
import geneticAlgorithm
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
        self.formLayout.addRow(self.numGenLabel, self.generationSpin)
        self.formLayout.addRow(self.tournamentSizeLabel, tournament_field)
        self.formLayout.addRow(self.elitismSizeLabel, self.elitismSpin)
        self.formLayout.addRow(self.mutationProbabilityLabel, self.mutationLineEdit)

        self.crossoverSmilesRadio = QRadioButton("SMILES", application)
        self.crossoverSelfiesRadio = QRadioButton("SELFIES", application)
        self.crossoverSelfiesRadio.setChecked(True)
        self.crossoverGroup = QButtonGroup(application)
        self.crossoverGroup.addButton(self.crossoverSmilesRadio, 0)
        self.crossoverGroup.addButton(self.crossoverSelfiesRadio, 1)

        crossover_field = QWidget(application)
        crossover_row = QHBoxLayout(crossover_field)
        crossover_row.setContentsMargins(0, 0, 0, 0)
        crossover_row.setSpacing(16)
        crossover_row.addWidget(self.crossoverSmilesRadio, 0, Qt.AlignLeft)
        crossover_row.addWidget(self.crossoverSelfiesRadio, 0, Qt.AlignLeft)
        crossover_row.addStretch(1)
        self.formLayout.addRow(self.crossoverLabel, crossover_field)

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
        self.container.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

    def onRouletteChanged(self, state):
        if state == 2:
            self.tournamentSpin.setDisabled(True)
        else:
            self.tournamentSpin.setDisabled(False)

    def getGAParametersWidget(self):
        return self.container

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
        self.application.gaConfig.useSelfiesCrossover = self.crossoverSelfiesRadio.isChecked()

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

        moleculeBoxes.rightCont3.setLayout(moleculeBoxes.evolutionControlsLayout)
        moleculeBoxes.rightCont3.setMinimumWidth(380)
        moleculeBoxes.rightCont3.setMaximumWidth(520)
        moleculeBoxes.rightCont3.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

        # Show evolution / progress on stage 3 before the GA blocks the event loop.
        self.application.show_stage_3()
        self.application.viewEvolutionStageButton.setVisible(True)

        self.rouletteCheckBox.setDisabled(True)
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
                border: 2px solid gray;
            }
            QCheckBox::indicator:unchecked {
                background-color: lightgray;
                border: 2px solid gray;
            }
        """)
        self.generationSpin.setDisabled(True)
        self.tournamentSpin.setDisabled(True)
        self.elitismSpin.setDisabled(True)
        self.mutationLineEdit.setDisabled(True)
        self.crossoverSmilesRadio.setDisabled(True)
        self.crossoverSelfiesRadio.setDisabled(True)

        self.application.sbmtBtn.setDisabled(True)
        self.application.resBtn.setDisabled(True)

        self.application.blockTransfer = True
        self.application._cancel_evolution = False

        # Paint progress widgets once before the long-running GA blocks the event loop.
        QApplication.processEvents()

        moleculeBoxes.newGenerationMolecules = geneticAlgorithm.geneticAlgorithm(
            moleculeBoxes.selectedMolecules,
            True,
            self.application.gaConfig.generations,
            self.application.gaConfig.rouletteSelection,
            self.application.gaConfig.tournamentSize,
            self.application.gaConfig.elitismSize,
            self.application.gaConfig.mutationProbability,
            self.application.gaConfig.useSelfiesCrossover,
            self.application.getMutationInfo(),
            moleculeBoxes.individualLabel,
            moleculeBoxes.individualProgress,
            cancel_check=lambda: getattr(self.application, "_cancel_evolution", False),
            on_generation_start=moleculeBoxes._live_generation_on_start,
            on_new_individual=moleculeBoxes._live_generation_on_new_individual,
        )

        if getattr(self.application, "_cancel_evolution", False):
            return

        moleculeBoxes.loadNewGeneration(tuple(self.application.sliderValues))
        
        

