from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QLabel, QSpinBox, QWidget, QPushButton, QLineEdit, QCheckBox
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon
import geneticAlgorithm

class GAParameters:
    def __init__(self, application):
        self.application = application
        self.mainLayout = QHBoxLayout()
        self.leftLayout = QVBoxLayout()
        self.gaLabel = QLabel("Parameters for genetic algorithm:", application)
        self.gaLabel.setStyleSheet("font-size: 17px; font-weight: bold; border: none; margin-bottom: 10px;")
        self.numGenLabel = QLabel("Number of generations: ", application)
        self.tournamentSizeLabel = QLabel("Tournament size: ", application)
        self.elitismSizeLabel = QLabel("Elitism size: ", application)
        self.mutationProbabilityLabel = QLabel("Mutation probability: ", application)
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
                background-color: #45a049; 
                border: 2px solid green;
            }
            QCheckBox::indicator:unchecked {
                background-color: white;
                border: 2px solid #777;
            }
        """)

        self.rouletteCheckBox.stateChanged.connect(self.onRouletteChanged)

        self.tournamentSpin = QSpinBox(application)
        self.tournamentSpin.setMaximum(100)
        self.tournamentSpin.setFixedWidth(70)
        self.tournamentSpin.setValue(7)
        self.elitismSpin = QSpinBox(application)
        self.elitismSpin.setMaximum(100)
        self.elitismSpin.setFixedWidth(70)
        self.elitismSpin.setValue(10)
        self.mutationLineEdit = QLineEdit(application)
        self.mutationLineEdit.setFixedWidth(70)
        self.mutationLineEdit.setText("0.05")

        self.h1 = QHBoxLayout()
        self.h1.addWidget(self.numGenLabel)
        self.h1.addWidget(self.generationSpin)
        self.h1Cont = QWidget()
        self.h1Cont.setLayout(self.h1)
        self.h1Cont.setFixedSize(300, 37)

        self.h2 = QHBoxLayout()
        self.h2.addWidget(self.tournamentSizeLabel)
        self.h2.addWidget(self.tournamentSpin)
        self.h2.addSpacing(20)
        self.h2.addWidget(self.rouletteCheckBox)
        self.h2Cont = QWidget()
        self.h2Cont.setLayout(self.h2)
        self.h2Cont.setFixedSize(540, 37)

        self.h3 = QHBoxLayout()
        self.h3.addWidget(self.elitismSizeLabel)
        self.h3.addWidget(self.elitismSpin)
        self.h3Cont = QWidget()
        self.h3Cont.setLayout(self.h3)
        self.h3Cont.setFixedSize(300, 37)

        self.h4 = QHBoxLayout()
        self.h4.addWidget(self.mutationProbabilityLabel)
        self.h4.addWidget(self.mutationLineEdit)
        self.h4Cont = QWidget()
        self.h4Cont.setLayout(self.h4)
        self.h4Cont.setFixedSize(300, 37)

        self.leftLayout.addWidget(self.gaLabel)
        self.leftLayout.addWidget(self.h1Cont)
        self.leftLayout.addWidget(self.h2Cont)
        self.leftLayout.addWidget(self.h3Cont)
        self.leftLayout.addWidget(self.h4Cont)
        self.leftLayout.setAlignment(Qt.AlignLeft)

        self.rightLayout = QVBoxLayout()
        self.rightLayout.setAlignment(Qt.AlignCenter)

        self.launchButton = QPushButton("Launch search!", application)
        self.launchButton.setFixedWidth(200)
        # self.launchButton.setIcon(QIcon('icon.png'))
        # self.launchButton.setIconSize(QSize(32, 32))
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

        self.rightLayout.addWidget(self.launchButton)

        self.leftWrapper = QWidget()
        self.leftWrapper.setLayout(self.leftLayout)

        self.rightWrapper = QWidget()
        self.rightWrapper.setLayout(self.rightLayout)
        self.rightWrapper.setFixedWidth(230)

        self.mainLayout.addWidget(self.leftWrapper)
        self.mainLayout.addWidget(self.rightWrapper)
        self.container = QWidget()
        self.container.setLayout(self.mainLayout)
        self.container.setFixedSize(760, 200)

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
        moleculeBoxes.generateButton.setDisabled(False)
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
        moleculeBoxes.finalButton.setDisabled(False)
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
        moleculeBoxes.saveButton.setDisabled(False)
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
        rouletteSelection = self.rouletteCheckBox.isChecked()
        numberOfGenerations = self.generationSpin.value()
        tournamentSize = self.tournamentSpin.value()
        elitismSize = self.elitismSpin.value()
        mutationProbability = float(self.mutationLineEdit.text())

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

        self.application.blockTransfer = True

        geneticAlgorithm.geneticAlgorithm(
            moleculeBoxes.selectedMolecules,
            True,
            numberOfGenerations,
            rouletteSelection,
            tournamentSize,
            elitismSize,
            mutationProbability
        )
        
        

