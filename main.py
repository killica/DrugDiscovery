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

class Application(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Drug Discovery')
        self.resize(800, 600)

        # self.selectionLabel = QLabel("Select molecules for the first generation:", self)

        #self.selectionLabel.setStyleSheet("font-size: 20px; font-weight: bold; padding-bottom: 10px;")

        self.mainLayout = QHBoxLayout()
        self.leftLayout = QVBoxLayout()

        # self.scrollArea = QScrollArea()
        # self.scrollArea.setWidgetResizable(True)
        # self.scrollArea.setFixedSize(760, 290)
        # self.scrollWidget = QWidget()
        # self.gridLayout = QGridLayout()

        self.molecules = self.readMolecules()

        self.moleculeBoxes = MoleculeBoxes(self.molecules, self.width())
        self.newMoleculeForm = NewMoleculeForm(self)
        self.hyperParamLayout = HyperParameters(self)
        self.gaParameters = GAParameters(self)

        # h1 will contain form on the left and slide bars on the right
        self.cnt = QWidget()
        self.h1 = QHBoxLayout()
        self.h1.setSizeConstraint(760)

        self.h1.addWidget(self.newMoleculeForm.getForm())
        self.h1.addWidget(self.hyperParamLayout.getSlidersWidget())

        self.cnt.setLayout(self.h1)
        self.cnt.setFixedWidth(765)
        self.cnt.setFixedHeight(300)

        # self.scrollWidget.setLayout(self.gridLayout)
        # self.scrollArea.setWidget(self.scrollWidget)
        self.leftLayout.addWidget(self.moleculeBoxes.getSelectionWidget())
        self.leftLayout.addSpacing(30)
        self.leftLayout.addWidget(self.cnt)
        self.leftLayout.addSpacing(30)
        self.leftLayout.addWidget(self.gaParameters.getGAParametersWidget())
        

        self.leftWrapper = QWidget()
        self.leftWrapper.setLayout(self.leftLayout)
        self.leftWrapper.setFixedHeight(880)
        self.mainLayout.addWidget(self.leftWrapper)
        self.mainLayout.addWidget(QLabel("123123123", self))
        self.mainLayout.setAlignment(Qt.AlignTop)

        self.setLayout(self.mainLayout)
        self.setFixedSize(1700, 900)

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