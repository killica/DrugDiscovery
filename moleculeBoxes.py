import json
from PyQt5.QtWidgets import QGroupBox, QVBoxLayout, QLabel, QGraphicsDropShadowEffect
from PyQt5.QtGui import QImage, QPixmap, QColor
from rdkit import Chem
from rdkit.Chem import Draw
from fitness import Fitness

class MoleculeBoxes:
    def __init__(self, molecules, layout, windowWidth, scrollArea):
        self.molecules = molecules
        self.layout = layout
        self.windowWidth = windowWidth
        self.scrollArea = scrollArea
        self.boxWidth = 220
        self.columnsPerRow = self.windowWidth // self.boxWidth
        self.columnsPerRow = max(1, self.columnsPerRow)
        for index, (smiles, description, qed) in enumerate(self.molecules):
            row = index // self.columnsPerRow 
            col = index % self.columnsPerRow
            # fit = Fitness(Chem.MolFromSmiles(smiles))
            # qed = fit.qed()
            description += "\n" + "QED: " + str(round(qed, 4))
            moleculeBox = self.createMoleculeBox(smiles, description, qed)
            self.layout.addWidget(moleculeBox, row, col)

    def createMoleculeBox(self, smiles, description, qed):
        box = QGroupBox()
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
        description += "\n" + "QED: " + str(round(qed, 4))
        moleculeBox = self.createMoleculeBox(smiles, description, qed)
        self.layout.addWidget(moleculeBox, row, col)
        self.scrollArea.verticalScrollBar().setValue(self.scrollArea.verticalScrollBar().maximum())
        with open('molecules.json', 'r') as file:
            data = json.load(file)
        data.append({
            "SMILES": smiles,
            "Description": description,
            "QED": qed
        })
        with open('molecules.json', 'w') as file:
            json.dump(data, file, indent = 4)
    
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
