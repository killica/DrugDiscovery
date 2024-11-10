import json
from PyQt5.QtWidgets import QGroupBox, QVBoxLayout, QHBoxLayout, QLabel, QGraphicsDropShadowEffect, QWidget, QGridLayout, QScrollArea
from PyQt5.QtGui import QImage, QPixmap, QColor, QFont
from PyQt5.QtCore import QEvent, Qt
from rdkit import Chem
from rdkit.Chem import Draw
from fitness import Fitness

class ClickableGroupBox(QGroupBox):
    def __init__(self, moleculeBoxes, index, ind, parent=None):
        super().__init__(parent)
        self.moleculeBoxes = moleculeBoxes
        self.index = index
        self.ind = ind
        
    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            # print("LUJKO! " + self.layout().itemAt(1).widget().text())
            self.moleculeBoxes.removeBoxes()
            if self.ind == 0:
                self.moleculeBoxes.selectedMolecules.append(self.moleculeBoxes.molecules[self.index])
                self.moleculeBoxes.molecules.pop(self.index)
            elif self.ind == 1:
                self.moleculeBoxes.molecules.append(self.moleculeBoxes.selectedMolecules[self.index])
                self.moleculeBoxes.selectedMolecules.pop(self.index)
            self.moleculeBoxes.loadBoxes()
        super().mousePressEvent(event)

    def enterEvent(self, event):
        self.layout().itemAt(1).widget().setStyleSheet("color: darkgreen; font-weight: bold;")

    def leaveEvent(self, event):
        self.layout().itemAt(1).widget().setStyleSheet("color: black; font-weight: normal;")

class MoleculeBoxes(QWidget):
    def __init__(self, molecules, windowWidth):
        self.molecules = molecules
        self.windowWidth = windowWidth
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

        self.precedentLayout = QGridLayout()

        self.rightVBox1 = QVBoxLayout()

        self.precedentScrollArea = QScrollArea()
        self.precedentScrollArea.setWidgetResizable(True)
        self.precedentScrollArea.setFixedSize(760, 290)

        self.precedentScrollWidget = QWidget()
        self.precedentScrollWidget.setLayout(self.precedentLayout)
        self.precedentScrollArea.setWidget(self.precedentScrollWidget)

        self.precedentLabel = QLabel("1. generation")
        self.precedentLabel.setStyleSheet("font-size: 20px; font-weight: bold;")

        self.rightVBox1.addWidget(self.precedentLabel)
        self.rightVBox1.addWidget(self.precedentScrollArea)
        self.rightVBox1.setAlignment(Qt.AlignTop)

        self.rightCont1 = QWidget()
        self.rightCont1.setLayout(self.rightVBox1)
        self.rightCont1.setFixedSize(800, 350)

        self.loadBoxes()
        
    def loadBoxes(self):
        self.boxes = []
        self.selectedBoxes = []
        self.boxWidth = 220
        self.columnsPerRow = self.windowWidth // self.boxWidth
        self.columnsPerRow = max(1, self.columnsPerRow)
        for index, (smiles, description, qed) in enumerate(self.molecules):
            row = index // self.columnsPerRow 
            col = index % self.columnsPerRow
            # fit = Fitness(Chem.MolFromSmiles(smiles))
            # qed = fit.qed()
            description += "\n" + "QED: " + str(round(qed, 4))
            self.moleculeBox = self.createMoleculeBox(smiles, description, qed, index, 0)
            self.boxes.append(self.moleculeBox)
            self.layout.addWidget(self.moleculeBox, row, col)

        for index, (smiles, description, qed) in enumerate(self.selectedMolecules):
            row = index // self.columnsPerRow 
            col = index % self.columnsPerRow
            # fit = Fitness(Chem.MolFromSmiles(smiles))
            # qed = fit.qed()
            description += "\n" + "QED: " + str(round(qed, 4))
            self.selectedMoleculeBox = self.createMoleculeBox(smiles, description, qed, index, 1)
            self.selectedBoxes.append(self.selectedMoleculeBox)
            self.precedentLayout.addWidget(self.selectedMoleculeBox, row, col)

    def removeBoxes(self):
        for box in self.boxes:
            self.layout.removeWidget(box)
            box.deleteLater()
        self.boxes = []
        for selectedBox in self.selectedBoxes:
            self.precedentLayout.removeWidget(selectedBox)
            selectedBox.deleteLater()
        self.selectedBoxes = []

    def getSelectionWidget(self):
        return self.container

    def getPrecedentScrollArea(self):
        return self.rightCont1

    def createMoleculeBox(self, smiles, description, qed, index, ind):
        box = ClickableGroupBox(self, index, ind)
        boxLayout = QVBoxLayout()
        mol = Chem.MolFromSmiles(smiles)
        molImage = Draw.MolToImage(mol, size=(200, 200))
        qimage = QImage(molImage.tobytes(), molImage.width, molImage.height, molImage.width * 3, QImage.Format_RGB888)
        pixmap = QPixmap(qimage)
        imageLabel = QLabel()
        imageLabel.setPixmap(pixmap)
        descriptionLabel = QLabel(description)
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
        fit = Fitness(Chem.MolFromSmiles(smiles))
        qed = fit.qed()
        
        self.removeBoxes()
        self.molecules.append((smiles, description, qed))
        with open('molecules.json', 'r') as file:
            data = json.load(file)
        data.append({
            "SMILES": smiles,
            "Description": description,
            "QED": qed
        })
        with open('molecules.json', 'w') as file:
            json.dump(data, file, indent = 4)
        self.loadBoxes()
        
        # self.addMoleculeBox(smiles, description, qed, ind)
