from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QLabel, QSlider, QWidget, QPushButton
from PyQt5.QtCore import Qt

class HyperParameters:
    def __init__(self, application):
        self.names = ['MW', 'ALOGP', 'HBA', 'HBD', 'PSA', 'ROTB', 'AROM', 'ALERTS']
        self.defaultValues = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95]
        self.hBoxes = []
        self.hyperParamLayout = QVBoxLayout()
        self.adjustLabel = QLabel("Adjust hyperparameter's weights", application)
        self.adjustLabel.setStyleSheet("font-size: 17px; font-weight: bold; border: none; margin-top: 5px;")

        self.resetButton = QPushButton("Reset to defaults", application)
        self.resetButton.clicked.connect(self.onResetButtonClicked)
        self.resetButton.setFixedWidth(120)
        self.resetButton.setStyleSheet("""
            QPushButton {
                background-color: #6495ED;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #0047AB;
            }
        """)
        self.hCont = QHBoxLayout()
        self.hCont.addWidget(self.adjustLabel)
        self.hCont.addWidget(self.resetButton)

        self.hyperParamLayout.addLayout(self.hCont)
        self.hyperParamLayout.addSpacing(20)
        for i, (name, default) in enumerate(zip(self.names, self.defaultValues)):
            self.hBox = QHBoxLayout()

            self.paramLabel = QLabel(name + ": " + str(default), application)

            self.paramSlider = QSlider(Qt.Horizontal, application)
            self.paramSlider.setRange(0, 100)
            self.paramSlider.setTickInterval(1) 
            self.paramSlider.setFixedWidth(350)
            self.paramSlider.setValue(default * 100)
            self.paramSlider.valueChanged.connect(lambda value, idx = i: self.updateParamLabel(idx))
            self.paramSlider.setStyleSheet("""
                border: none;               
            """)

            self.hBox.addWidget(self.paramLabel)
            self.hBox.addWidget(self.paramSlider)
            self.hBoxes.append(self.hBox)
            self.hyperParamLayout.addLayout(self.hBox)

        self.container = QWidget()
        self.container.setLayout(self.hyperParamLayout)
        self.container.setStyleSheet("""  
            padding-top: 0px;            
        """)
            
    def updateParamLabel(self, idx):
        value = self.hBoxes[idx].itemAt(1).widget().value()
        name = self.names[idx]
        label = self.hBoxes[idx].itemAt(0).widget()
        label.setText(f"{name}: {value/100}")

    def onResetButtonClicked(self):
        for i in range(len(self.names)):
            self.hBoxes[i].itemAt(1).widget().setValue(self.defaultValues[i] * 100)



    def getSlidersWidget(self):
        return self.container