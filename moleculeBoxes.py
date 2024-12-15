import json
import copy
import re
from PyQt5.QtWidgets import QGroupBox, QVBoxLayout, QHBoxLayout, QLabel, QGraphicsDropShadowEffect, QWidget, QGridLayout, QScrollArea, QPushButton
from PyQt5.QtGui import QImage, QPixmap, QColor, QFont
from PyQt5.QtCore import QEvent, Qt
from rdkit import Chem
from rdkit.Chem import Draw
from fitness import Fitness
from individual import Individual
import geneticAlgorithm
from mutationInfo import MutationInfo
from datetime import datetime

class ClickableGroupBox(QGroupBox):
    def __init__(self, moleculeBoxes, index, ind, parent=None):
        super().__init__(parent)
        self.moleculeBoxes = moleculeBoxes
        self.index = index
        self.ind = ind
        
    def mousePressEvent(self, event):
        if self.moleculeBoxes.application.blockTransfer:
            return
        if self.ind != 0 and self.ind != 1:
            return
        if event.button() == Qt.LeftButton:
            self.moleculeBoxes.removeBoxes()
            self.moleculeBoxes.removeSelectedBoxes()
            if self.ind == 0:
                self.moleculeBoxes.selectedMolecules.append(self.moleculeBoxes.molecules[self.index])
                self.moleculeBoxes.molecules.pop(self.index)
            elif self.ind == 1:
                self.moleculeBoxes.molecules.append(self.moleculeBoxes.selectedMolecules[self.index])
                self.moleculeBoxes.selectedMolecules.pop(self.index)
            self.moleculeBoxes.loadBoxes(tuple(self.moleculeBoxes.application.sliderValues))
            self.moleculeBoxes.loadSelectedBoxes(tuple(self.moleculeBoxes.application.sliderValues))

        super().mousePressEvent(event)

    def enterEvent(self, event):
        if not self.moleculeBoxes.application.blockTransfer:
            self.layout().itemAt(1).widget().setStyleSheet("color: darkgreen; font-weight: bold;")

    def leaveEvent(self, event):
        if not self.moleculeBoxes.application.blockTransfer:
            self.layout().itemAt(1).widget().setStyleSheet("color: black; font-weight: normal;")

class MoleculeBoxes(QWidget):
    def __init__(self, application):
        super().__init__(application)
        self.application = application
        self.molecules = application.molecules
        self.windowWidth = application.width()
        self.boxWidth = 220
        self.columnsPerRow = 3
        self.layout = QGridLayout()

        self.vbox = QVBoxLayout()

        self.selectionLabel = QLabel("Select molecules for the first generation:")
        self.selectionLabel.setStyleSheet("font-size: 20px; font-weight: bold; padding-bottom: 10px;")

        self.scrollArea = QScrollArea()
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setFixedSize(760, 290)

        self.scrollWidget = QWidget()
        self.scrollWidget.setLayout(self.layout)
        self.scrollArea.setWidget(self.scrollWidget)

        self.vbox.addWidget(self.selectionLabel)
        self.vbox.addWidget(self.scrollArea)

        self.container = QWidget()
        self.container.setLayout(self.vbox)
        self.container.setFixedSize(760, 350)

        self.selectedMolecules = []
        self.newGenerationMolecules = []

        self.precedentLayout = QGridLayout()

        self.rightHBox1 = QHBoxLayout()

        self.precedentScrollArea = QScrollArea()
        self.precedentScrollArea.setWidgetResizable(True)
        self.precedentScrollArea.setFixedSize(760, 290)

        self.precedentScrollWidget = QWidget()
        self.precedentScrollWidget.setLayout(self.precedentLayout)
        self.precedentScrollArea.setWidget(self.precedentScrollWidget)

        self.precedentLabel = QLabel("1. generation")
        self.precedentLabel.setStyleSheet("font-size: 20px; font-weight: bold; color: darkgreen;")

        self.rightHBox1.addWidget(self.precedentScrollArea)
        self.rightHBox1.addSpacing(20)
        self.rightHBox1.addWidget(self.precedentLabel)

        self.rightHBox1.setAlignment(Qt.AlignTop)

        self.rightCont1 = QWidget()
        self.rightCont1.setLayout(self.rightHBox1)
        self.rightCont1.setFixedSize(920, 325)
        
        self.secondLayout = QGridLayout()

        self.rightHB = QHBoxLayout()

        self.secondScrollArea = QScrollArea()
        self.secondScrollArea.setWidgetResizable(True)
        self.secondScrollArea.setFixedSize(760, 290)

        self.secondScrollWidget = QWidget()
        self.secondScrollWidget.setLayout(self.secondLayout)
        self.secondScrollArea.setWidget(self.secondScrollWidget)

        self.secondLabel = QLabel("2. generation")
        self.secondLabel.setStyleSheet("font-size: 20px; font-weight: bold; color: darkgreen;")

        self.rightHB.addWidget(self.secondScrollArea)
        self.rightHB.addSpacing(20)
        self.rightHB.addWidget(self.secondLabel)

        self.rightHB.setAlignment(Qt.AlignTop)

        self.rightCont2 = QWidget()
        self.rightCont2.setLayout(self.rightHB)
        self.rightCont2.setFixedSize(920, 325)

        self.rightHBox2 = QHBoxLayout()

        self.rightVBox3 = QVBoxLayout()

        self.generateButton = QPushButton("Generate next", application)
        self.generateButton.setDisabled(True)
        self.generateButton.setFixedWidth(150)
        self.generateButton.clicked.connect(self.onGenerateButtonClicked)

        self.finalButton = QPushButton("Jump to final")
        self.finalButton.setDisabled(True)
        self.finalButton.setFixedWidth(150)
        self.finalButton.clicked.connect(self.onFinalButtonClicked)

        self.saveButton = QPushButton("Save the best")
        self.saveButton.setDisabled(True)
        self.saveButton.setFixedWidth(150)
        self.saveButton.clicked.connect(self.onSaveButtonClicked)

        self.saveLabel = QLabel("Saved!")
        self.saveLabel.setStyleSheet("color: transparent; font-style: italic;")

        self.saveBox = QHBoxLayout()
        self.saveBox.addSpacing(-7)
        self.saveBox.addWidget(self.saveButton)
        self.saveBox.addSpacing(10)
        self.saveBox.addWidget(self.saveLabel)
        # self.saveBox.setAlignment(Qt.AlignLeft)

        self.saveCnt = QWidget()
        self.saveCnt.setLayout(self.saveBox)
        self.saveCnt.setFixedWidth(350)

        self.restartButton = QPushButton("Restart analysis")
        self.restartButton.setDisabled(True)
        self.restartButton.setFixedWidth(150)
        self.restartButton.clicked.connect(self.onRestartButtonClicked)

        self.rightVBox3.addWidget(self.generateButton)
        self.rightVBox3.addSpacing(13)
        self.rightVBox3.addWidget(self.finalButton)
        self.rightVBox3.addWidget(self.saveCnt)
        self.rightVBox3.addWidget(self.restartButton)

        self.rightBtnCnt = QWidget()
        self.rightBtnCnt.setLayout(self.rightVBox3)
        self.rightBtnCnt.setFixedSize(350, 215)

        # self.bestBox = self.createMoleculeBox("CN(C)CCCN1C2=CC=CC=C2SC3=C1C=C(C=C3)Cl", "Chlorpromazine", 0.55, 0, -1)
        self.bestBox = self.createMoleculeBox("", "To be determined", 0.0, 0, -1)
        self.bestBox.setAlignment(Qt.AlignCenter)

        self.rightHBox2.addWidget(self.rightBtnCnt)
        self.rightHBox2.addWidget(self.bestBox)
        self.rightHBox2.setAlignment(Qt.AlignHCenter)

        self.rightCont3 = QWidget()
        self.rightCont3.setLayout(self.rightHBox2)
        self.rightCont3.setFixedSize(520, 270)

        self.loadBoxes()
        self.loadSelectedBoxes()
        
    def loadBoxes(self, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self.boxes = []
        self.molecules.sort(reverse=True)
        for index, individual in enumerate(self.molecules):
            row = index // 3
            col = index % 3
            smiles = individual.getSmiles()
            description = individual.getDescription()
            individual.setWeights(weights)
            qed = individual.getQED()
            self.moleculeBox = self.createMoleculeBox(smiles, description, qed, index, 0)
            self.boxes.append(self.moleculeBox)
            self.layout.addWidget(self.moleculeBox, row, col)

    def loadSelectedBoxes(self, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self.selectedBoxes = []
        self.selectedMolecules.sort(reverse=True)
        for index, individual in enumerate(self.selectedMolecules):
            row = index // 3
            col = index % 3
            smiles = individual.getSmiles()
            description = individual.getDescription()
            individual.setWeights(weights)
            qed = individual.getQED()
            self.selectedMoleculeBox = self.createMoleculeBox(smiles, description, qed, index, 1)
            self.selectedBoxes.append(self.selectedMoleculeBox)
            self.precedentLayout.addWidget(self.selectedMoleculeBox, row, col)

        # Move to bottom when a new moleecule box is added
        # self.scrollArea.verticalScrollBar().setValue(self.scrollArea.verticalScrollBar().maximum())
        # self.precedentScrollArea.verticalScrollBar().setValue(self.precedentScrollArea.verticalScrollBar().maximum())

    def loadNewGeneration(self, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self.newGenerationBoxes = []
        self.newGenerationMolecules.sort(reverse=True)
        for index, individual in enumerate(self.newGenerationMolecules):
            row = index // self.columnsPerRow 
            col = index % self.columnsPerRow
            smiles = individual.getSmiles()
            description = individual.getDescription()
            individual.setWeights(weights)
            qed = individual.getQED()
            self.newGenerationMoleculeBox = self.createMoleculeBox(smiles, description, qed, index, -1)
            self.newGenerationBoxes.append(self.newGenerationMoleculeBox)
            self.secondLayout.addWidget(self.newGenerationMoleculeBox, row, col)
        
        self.bestBox.deleteLater()
        self.bestBox = self.createMoleculeBox(self.newGenerationMolecules[0].getSmiles(), "Current best", self.newGenerationMolecules[0].getQED(), 0, -1)
        self.bestBox.setAlignment(Qt.AlignCenter)
        self.rightHBox2.addWidget(self.bestBox)

    def removeBoxes(self):
        for box in self.boxes:
            self.layout.removeWidget(box)
            box.deleteLater()
        self.boxes = []

    def removeSelectedBoxes(self):
        for selectedBox in self.selectedBoxes:
            self.precedentLayout.removeWidget(selectedBox)
            selectedBox.deleteLater()
        self.selectedBoxes = []

    def removeNewGenerationBoxes(self):
        for newGenerationBox in self.newGenerationBoxes:
            self.secondLayout.removeWidget(newGenerationBox)
            newGenerationBox.deleteLater()
        self.newGenerationBoxes = []

    def getSelectionWidget(self):
        return self.container

    def getPrecedentScrollArea(self):
        return self.rightCont1

    def getSecondScrollArea(self):
        return self.rightCont2

    def getBest(self):
        return self.rightCont3

    def createMoleculeBox(self, smiles, description, qed, index, ind):
        box = ClickableGroupBox(self, index, ind)
        box.setFixedWidth(230)
        boxLayout = QVBoxLayout()
        mol = Chem.MolFromSmiles(smiles)
        molImage = Draw.MolToImage(mol, size=(200, 200))
        qimage = QImage(molImage.tobytes(), molImage.width, molImage.height, molImage.width * 3, QImage.Format_RGB888)
        pixmap = QPixmap(qimage)
        imageLabel = QLabel()
        imageLabel.setPixmap(pixmap)
        descriptionLabel = QLabel(description + "\n" + "QED: " + str(round(qed, 4)))
        boxLayout.addWidget(imageLabel)
        boxLayout.addWidget(descriptionLabel)
        box.setLayout(boxLayout)

        # Formula to determine the shade of green - darker the shade, greater the QED
        # box.setStyleSheet(f"QGroupBox {{ background-color: rgb({0}, {int(255-(qed*155))}, {0}); border: 2px solid green; padding: 10px; }}")
        
        box.setStyleSheet(f"""
            QGroupBox {{
                background-color: rgb(255, 255, 255);   
                border: 2px solid green;                
                border-radius: 15px;                    
                padding: 10px;                          
            }}
        """)

        shadowEffect = QGraphicsDropShadowEffect()
        shadowEffect.setOffset(5, 5)              
        shadowEffect.setBlurRadius(15)            
        shadowEffect.setColor(QColor(0, 0, 0, 160))
        box.setGraphicsEffect(shadowEffect)
        return box
    
    def addToCatalogue(self, smiles, description):
        molecule = None
        try:
            molecule = Chem.MolFromSmiles(smiles)
        except ValueError:
            return
        if molecule is None:
            return
        if not description:
            description = "Unknown"
        
        self.removeBoxes()
        self.removeSelectedBoxes()
        self.molecules.append(Individual(smiles, description, self.application.sliderValues))
        with open('molecules.json', 'r') as file:
            data = json.load(file)
        data.append({
            "SMILES": smiles,
            "Description": description
        })
        with open('molecules.json', 'w') as file:
            json.dump(data, file, indent = 4)
        self.loadBoxes()
        self.loadSelectedBoxes()

    def onGenerateButtonClicked(self):
        self.saveLabel.setStyleSheet("color: transparent; font-style: italic;")
        self.removeSelectedBoxes()
        self.selectedMolecules = []
        for ind in self.newGenerationMolecules:
            self.selectedMolecules.append(Individual(ind.getSmiles(), ind.getDescription()))
        self.loadSelectedBoxes(tuple(self.application.sliderValues))
        self.removeNewGenerationBoxes()

        self.newGenerationMolecules = geneticAlgorithm.geneticAlgorithm(
            self.selectedMolecules,
            True,
            self.application.numberOfGenerations,
            self.application.rouletteSelection,
            self.application.tournamentSize,
            self.application.elitismSize,
            self.application.mutationProbability,
            self.application.mi
        )

        self.loadNewGeneration(tuple(self.application.sliderValues))
        labelText = self.secondLabel.text()
        self.precedentLabel.setText(labelText)
        # Regular expression to match a number at the start of the string
        match = re.match(r'^\d+', labelText)
        # Check if a match was found and extract the number
        labelNumber = int(match.group(0)) + 1
        self.secondLabel.setText(str(labelNumber) + ". generation")

        
    def onFinalButtonClicked(self):
        self.saveLabel.setStyleSheet("color: transparent; font-style: italic;")
        labelText = self.secondLabel.text()
        # Regular expression to match a number at the start of the string
        match = re.match(r'^\d+', labelText)
        # Check if a match was found and extract the number
        labelNumber = int(match.group(0))
        for _ in range(self.application.numberOfGenerations - labelNumber):
            self.onGenerateButtonClicked()


    def onSaveButtonClicked(self):
        formattedDatetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open('best_candidate_molecules.txt', 'a') as candidatesFile:
            candidatesFile.write(f"SMILES: {self.newGenerationMolecules[0].getSmiles()}\nQED: {round(self.newGenerationMolecules[0].getQED(), 4)}\nDate created: {formattedDatetime}\n-------------------------------------------\n")
        self.saveLabel.setStyleSheet("color: green; font-style: italic;")

    def onRestartButtonClicked(self):
        self.saveLabel.setStyleSheet("color: transparent; font-style: italic;")

