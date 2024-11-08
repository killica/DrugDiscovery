import sys
from PyQt5.QtWidgets import QApplication, QWidget, QDesktopWidget, QPushButton, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QLineEdit, QScrollArea, QGroupBox
from PyQt5.QtGui import QImage, QPixmap
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw
from moleculeBoxes import MoleculeBoxes
from insertMolecule import NewMoleculeForm

class Application(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Drug Discovery')
        self.resize(800, 600)
        self.center()

        self.selectionLabel = QLabel("Select molecules for the first generation:", self)

        self.selectionLabel.setStyleSheet("font-size: 20px; font-weight: bold; padding-bottom: 10px;")

        mainLayout = QVBoxLayout()

        scrollArea = QScrollArea()
        scrollWidget = QWidget()
        self.gridLayout = QGridLayout()

        self.molecules = [
            ('CC(=O)Oc1ccccc1C(=O)O', 'Aspirin (acetylsalicylic acid)'),
            ('CC(=O)N1CCCC1C(=O)O', 'Ibuprofen'),
            ('C1=CC=C(C=C1)C(=O)O', 'Benzoic Acid'),
            ('CC(=O)OCC(C(=O)O)C', 'Malonic Acid'),
            ('CC1=CC=CC=C1C(=O)O', 'Salicylic Acid'),
            ('CC(=O)Nc1ccc(C)cc1', 'Acetanilide'),
            ('C', 'Methane')
        ]

        self.moleculeBoxes = MoleculeBoxes(self.molecules, self.gridLayout, self.width())
        self.newMoleculeForm = NewMoleculeForm(self)

        scrollWidget.setLayout(self.gridLayout)
        scrollArea.setWidget(scrollWidget)
        mainLayout.addWidget(self.selectionLabel)
        mainLayout.addWidget(scrollArea)
        mainLayout.addLayout(self.newMoleculeForm.getLayout())
        mainLayout.addWidget(self.newMoleculeForm.getSubmitButton())

        self.setLayout(mainLayout)
        self.show()

    def onSubmitButtonClicked(self):
        smiles = self.newMoleculeForm.getInputSmilesText()
        description = self.newMoleculeForm.getInputDescriptionText()
        self.moleculeBoxes.addToCatalogue(smiles, description)

    def center(self):
        screen = QDesktopWidget().screenGeometry()
        windowSize = self.geometry()
        self.move(
            (screen.width() - windowSize.width()) // 2,
            (screen.height() - windowSize.height()) // 2
        )

    def resizeEvent(self, event):
        self.clearLayout(self.gridLayout)
        self.moleculeBoxes = MoleculeBoxes(self.molecules, self.gridLayout, self.width())
        super().resizeEvent(event)

    def clearLayout(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                self.clearLayout(item.layout())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = Application()
    sys.exit(app.exec_())