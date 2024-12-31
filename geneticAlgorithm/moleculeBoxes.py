import json
import copy
import re
from PyQt5.QtWidgets import QGroupBox, QVBoxLayout, QHBoxLayout, QLabel, QMessageBox, QGraphicsDropShadowEffect, QWidget, QGridLayout, QScrollArea, QPushButton, QProgressBar
from PyQt5.QtGui import QImage, QPixmap, QColor, QFont
from PyQt5.QtCore import QEvent, Qt
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors
from rdkit.DataStructs import FingerprintSimilarity
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
        self.precedentLabel.setFixedWidth(140)
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
        self.secondLabel.setFixedWidth(140)
        self.secondLabel.setStyleSheet("font-size: 20px; font-weight: bold; color: darkgreen;")

        self.rightHB.addWidget(self.secondScrollArea)
        self.rightHB.addSpacing(20)
        self.rightHB.addWidget(self.secondLabel)

        self.rightHB.setAlignment(Qt.AlignTop)

        self.rightCont2 = QWidget()
        self.rightCont2.setLayout(self.rightHB)
        self.rightCont2.setFixedSize(920, 325)

        self.rightHBox2 = QHBoxLayout()
        self.rightCont3 = QWidget()

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

    def loadNewGeneration(self, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self.newGenerationBoxes = []
        if len(self.newGenerationMolecules) > 0:
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
                self.rightHBox2.insertWidget(1, self.bestBox)
        # else:
        #     self.bestBox.deleteLater()
        #     self.bestBox = self.createMoleculeBox("", "To be determined", 0.0, 0, -1)
        #     self.bestBox.setAlignment(Qt.AlignCenter)
        #     self.rightHBox2.insertWidget(1, self.bestBox)
       

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
            self.application.mi,
            self.individualLabel,
            self.individualProgress
        )

        self.loadNewGeneration(tuple(self.application.sliderValues))
        labelText = self.secondLabel.text()
        self.precedentLabel.setText(labelText)
        # Regular expression to match a number at the start of the string
        match = re.match(r'^\d+', labelText)
        # Check if a match was found and extract the number
        labelNumber = int(match.group(0)) + 1
        self.secondLabel.setText(str(labelNumber) + ". generation")
        self.generationLabel.setText(f"Generation: {labelNumber}/{self.application.numberOfGenerations}")
        self.generationProgress.setValue(labelNumber)

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
        with open('results/best_candidate_molecules.txt', 'a') as candidatesFile:
            candidatesFile.write(f"SMILES: {self.newGenerationMolecules[0].getSmiles()}\nQED: {round(self.newGenerationMolecules[0].getQED(), 4)}\nDate created: {formattedDatetime}\n-------------------------------------------\n")
       
        self.saveLabel.setStyleSheet("color: green; font-style: italic;")

        # Create a QMessageBox
        msg_box = QMessageBox(self)
        # Set the icon for the dialog
        msg_box.setIcon(QMessageBox.Question)
        # Set the window title
        msg_box.setWindowTitle("Population diversity")
        # Set the message in the dialog
        msg_box.setText("Do you want to calculate Tanimoto similarity coefficient for current generation?")
        # Add Yes and No buttons
        msg_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        # Show the message box and capture the response
        response = msg_box.exec_()
        if response == QMessageBox.Yes:
            self.tanimoto()
       

    def tanimoto(self):
        smilesList = [s.getSmiles() for s in self.newGenerationMolecules]
        mols = [Chem.MolFromSmiles(s) for s in smilesList]
        fingerprints = [rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits = 2048) for mol in mols]
        # Calculate pairwise Tanimoto similarity
        similarities = []
        for i in range(len(fingerprints)):
            for j in range(i+1, len(fingerprints)):
                similarity = FingerprintSimilarity(fingerprints[i], fingerprints[j])
                similarities.append(similarity)
        with open('results/tanimoto.txt', 'a') as tanimotoFile:
            tanimotoFile.write("[")
            for sim in similarities:
                tanimotoFile.write(f"{sim}, ")
            tanimotoFile.write("]\n------------------------------------\n")

    def onRestartButtonClicked(self):
        self.saveLabel.setStyleSheet("color: transparent; font-style: italic;")
        self.generationLabel.setText(f"Generation: 1/{self.application.numberOfGenerations}")
        self.generationProgress.setValue(1)
        self.individualLabel.setText(f"Individual: 0/{len(self.selectedMolecules)}")
        self.individualProgress.setValue(0)
        self.removeBoxes()
        self.removeSelectedBoxes()
        self.removeNewGenerationBoxes()
        self.molecules = self.application.readMolecules()
        self.loadBoxes(self.application.sliderValues)
        self.application.molecules = []
        self.selectedMolecules = []
        self.newGenerationMolecules = []
        self.bestBox.deleteLater()
        self.generateButton.setDisabled(True)
        self.generateButton.setStyleSheet("Color: #757575;")
        self.finalButton.setDisabled(True)
        self.finalButton.setStyleSheet("Color: #757575;")
        self.saveButton.setDisabled(True)
        self.saveButton.setStyleSheet("Color: #757575;")
        self.restartButton.setDisabled(True)
        self.restartButton.setStyleSheet("Color: #757575;")
        self.loadNewGeneration()
        self.precedentLabel.setText("1. generation")
        self.secondLabel.setText("2. generation")
        self.application.gaParameters.launchButton.setDisabled(False)
        self.application.gaParameters.launchButton.setStyleSheet("""
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
        self.application.gaParameters.rouletteCheckBox.setDisabled(False)
        self.application.gaParameters.rouletteCheckBox.setStyleSheet("""
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
                background-color: white;
                border: 2px solid gray;
            }
        """)
        self.application.gaParameters.generationSpin.setDisabled(False)
        self.application.gaParameters.tournamentSpin.setDisabled(False)
        self.application.gaParameters.elitismSpin.setDisabled(False)
        self.application.gaParameters.mutationLineEdit.setDisabled(False)
        self.application.sbmtBtn.setDisabled(False)
        self.application.resBtn.setDisabled(False)
        self.application.blockTransfer = False
        while self.rightHBox2.count():
            item = self.rightHBox2.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()