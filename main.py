import sys
import json
from PyQt5.QtWidgets import QApplication, QWidget, QSlider, QDesktopWidget, QPushButton, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QLineEdit, QScrollArea, QGroupBox
from PyQt5.QtGui import QImage, QPixmap
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw
from moleculeBoxes import MoleculeBoxes
from insertMolecule import NewMoleculeForm
from hyperParameters import HyperParameters

class Application(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Drug Discovery')
        self.resize(800, 600)
        self.center()

        self.selectionLabel = QLabel("Select molecules for the first generation:", self)

        self.selectionLabel.setStyleSheet("font-size: 20px; font-weight: bold; padding-bottom: 10px;")

        self.mainLayout = QVBoxLayout()

        self.scrollArea = QScrollArea()
        self.scrollArea.setWidgetResizable(False)
        self.scrollArea.setFixedSize(760, 290)
        self.scrollWidget = QWidget()
        self.gridLayout = QGridLayout()

        self.molecules = self.readMolecules()

        self.moleculeBoxes = MoleculeBoxes(self.molecules, self.gridLayout, self.width(), self.scrollArea)
        self.newMoleculeForm = NewMoleculeForm(self)
        self.hyperParamLayout = HyperParameters(self)


        # h1 will contain form on the left and slide bars on the right
        self.cnt = QWidget()
        self.h1 = QHBoxLayout()
        self.h1.setSizeConstraint(760)

        self.h1.addWidget(self.newMoleculeForm.getForm())
        self.h1.addWidget(self.hyperParamLayout.getSlidersWidget())

        self.cnt.setLayout(self.h1)
        self.cnt.setFixedWidth(765)

        self.scrollWidget.setLayout(self.gridLayout)
        self.scrollArea.setWidget(self.scrollWidget)
        self.mainLayout.addWidget(self.selectionLabel)
        self.mainLayout.addWidget(self.scrollArea)
        self.mainLayout.addSpacing(30)
        self.mainLayout.addWidget(self.cnt)

        self.setLayout(self.mainLayout)
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

    # def resizeEvent(self, event):
    #     self.clearLayout(self.gridLayout)
    #     self.moleculeBoxes = MoleculeBoxes(self.molecules, self.gridLayout, self.width(), self.scrollArea)
    #     super().resizeEvent(event)

    def clearLayout(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                self.clearLayout(item.layout())

    def readMolecules(self):
        with open('molecules.json', 'r') as file:
            data = json.load(file)

        return [(item['SMILES'], item['Description'], item['QED']) for item in data]


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = Application()
    sys.exit(app.exec_())