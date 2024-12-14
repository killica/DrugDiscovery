import sys
import json
from PyQt5.QtWidgets import QApplication, QWidget, QSlider, QDesktopWidget, QPushButton, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QLineEdit, QScrollArea, QGroupBox
from PyQt5.QtGui import QImage, QPixmap, QPainter, QColor, QPen
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw
from moleculeBoxes import MoleculeBoxes
from insertMolecule import NewMoleculeForm
from hyperParameters import HyperParameters
from gaParameters import GAParameters
from individual import Individual

class Application(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Drug Discovery')
        self.resize(800, 600)

        self.mainLayout = QHBoxLayout()
        self.leftLayout = QVBoxLayout()

        self.molecules = self.readMolecules()
        self.sliderValues = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95]

        # Allow transfering molecule boxes between scroll areas
        self.blockTransfer = False

        self.moleculeBoxes = MoleculeBoxes(self)
        self.newMoleculeForm = NewMoleculeForm(self)
        self.hyperParamLayout = HyperParameters(self)
        self.gaParameters = GAParameters(self)

        self.sbmtBtn = self.newMoleculeForm.submitButton
        self.resBtn = self.hyperParamLayout.resetButton

        self.cnt = QWidget()
        self.h1 = QHBoxLayout()
        self.h1.setSizeConstraint(760)

        self.h1.addWidget(self.newMoleculeForm.getForm())
        self.h1.addWidget(self.hyperParamLayout.getSlidersWidget())

        self.cnt.setLayout(self.h1)
        self.cnt.setFixedWidth(765)
        self.cnt.setFixedHeight(300)

        self.leftLayout.addWidget(self.moleculeBoxes.getSelectionWidget())
        self.leftLayout.addSpacing(70)
        self.leftLayout.addWidget(self.cnt)
        self.leftLayout.addSpacing(30)
        self.leftLayout.addWidget(self.gaParameters.getGAParametersWidget())
        
        self.leftWrapper = QWidget()
        self.leftWrapper.setLayout(self.leftLayout)
        self.leftWrapper.setFixedHeight(880)
        self.mainLayout.addWidget(self.leftWrapper)

        self.rightLayout = QVBoxLayout()
        self.rightLayout.addWidget(self.moleculeBoxes.getPrecedentScrollArea())
        self.rightLayout.addWidget(self.moleculeBoxes.getSecondScrollArea())
        self.rightLayout.addWidget(self.moleculeBoxes.getBest())

        self.rightWrapper = QWidget()
        self.rightWrapper.setLayout(self.rightLayout)
        self.rightWrapper.setFixedHeight(880)
        self.mainLayout.addWidget(self.rightWrapper)
        self.mainLayout.setAlignment(Qt.AlignTop)

        self.setLayout(self.mainLayout)
        self.setFixedSize(1750, 900)

        self.show()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setPen(QPen(QColor(128, 128, 128), 2))
        painter.drawLine(10, 685, 800, 685)
        painter.drawLine(800, 20, 800, 880)
        painter.end()

    def onSubmitButtonClicked(self):
        smiles = self.newMoleculeForm.getInputSmilesText()
        description = self.newMoleculeForm.getInputDescriptionText()
        self.moleculeBoxes.addToCatalogue(smiles, description)

    def readMolecules(self):
        with open('molecules.json', 'r') as file:
            data = json.load(file)
        return [Individual(item['SMILES'], item['Description']) for item in data]

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = Application()
    sys.exit(app.exec_())