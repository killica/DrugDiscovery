import json
from PyQt5.QtWidgets import QGroupBox, QVBoxLayout, QLabel, QGraphicsDropShadowEffect, QWidget, QGridLayout, QScrollArea
from PyQt5.QtGui import QImage, QPixmap, QColor, QFont
from PyQt5.QtCore import QEvent, Qt
from rdkit import Chem
from rdkit.Chem import Draw
from fitness import Fitness

class ClickableGroupBox(QGroupBox):
    def __init__(self, moleculeBoxes, index, parent=None):
        super().__init__(parent)
        self.moleculeBoxes = moleculeBoxes
        self.index = index
        
    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            # print("LUJKO! " + self.layout().itemAt(1).widget().text())
            self.moleculeBoxes.removeBoxes()
            self.moleculeBoxes.molecules.pop(self.index)
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

        self.loadBoxes()
        
    def loadBoxes(self):
        self.boxes = []
        self.boxWidth = 220
        self.columnsPerRow = self.windowWidth // self.boxWidth
        self.columnsPerRow = max(1, self.columnsPerRow)
        for index, (smiles, description, qed) in enumerate(self.molecules):
            row = index // self.columnsPerRow 
            col = index % self.columnsPerRow
            # fit = Fitness(Chem.MolFromSmiles(smiles))
            # qed = fit.qed()
            description += "\n" + "QED: " + str(round(qed, 4))
            self.moleculeBox = self.createMoleculeBox(smiles, description, qed, index)
            self.boxes.append(self.moleculeBox)
            self.layout.addWidget(self.moleculeBox, row, col)

    def removeBoxes(self):
        for box in self.boxes:
            self.layout.removeWidget(box)
            box.deleteLater()
        self.boxes = []

    def getSelectionWidget(self):
        return self.container

    def createMoleculeBox(self, smiles, description, qed, index):
        box = ClickableGroupBox(self, index)
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

    def addMoleculeBox(self, smiles, description):
        index = len(self.molecules) - 1
        row = index // self.columnsPerRow
        col = index % self.columnsPerRow
        fit = Fitness(Chem.MolFromSmiles(smiles))
        qed = fit.qed()
        with open('molecules.json', 'r') as file:
            data = json.load(file)
        data.append({
            "SMILES": smiles,
            "Description": description,
            "QED": qed
        })
        with open('molecules.json', 'w') as file:
            json.dump(data, file, indent = 4)
        description += "\n" + "QED: " + str(round(qed, 4))
        moleculeBox = self.createMoleculeBox(smiles, description, qed)
        self.layout.addWidget(moleculeBox, row, col)
        self.scrollArea.verticalScrollBar().setValue(self.scrollArea.verticalScrollBar().maximum())
    
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
        self.molecules.append((smiles, description))
        self.addMoleculeBox(smiles, description)
