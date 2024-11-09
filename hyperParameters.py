from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QLabel, QSlider, QWidget
from PyQt5.QtCore import Qt

class HyperParameters:
    def __init__(self, application):
        self.names = ['MW', 'ALOGP', 'HBA', 'HBD', 'PSA', 'ROTB', 'AROM', 'ALERTS']
        self.defaultValues = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95]
        self.hBoxes = []
        self.hyperParamLayout = QVBoxLayout()
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
            border-bottom: 2px solid gray;     
            padding-top: 0px;            
        """)
        #self.container.setFixedWidth(900)
            
    def updateParamLabel(self, idx):
        # paramLabel.setText(f"MW: {value/100}")
        value = self.hBoxes[idx].itemAt(1).widget().value()
        name = self.names[idx]
        label = self.hBoxes[idx].itemAt(0).widget()
        label.setText(f"{name}: {value/100}")


    def getSlidersWidget(self):
        return self.container