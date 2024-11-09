from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QLabel, QLineEdit, QWidget
from PyQt5.QtCore import Qt

class NewMoleculeForm:
    def __init__(self, application):
        self.customMoleculeTitle = QLabel("Or, add your own molecule:", application)
        # self.customMoleculeTitle.setAlignment(Qt.AlignCenter)

        self.customMoleculeTitle.setStyleSheet("font-size: 17px; font-weight: bold; border: none;")

        self.inputSMILES = QLineEdit(application)
        self.inputSMILES.setPlaceholderText("Enter SMILES format...")
        self.inputSMILES.setFixedWidth(250)
        self.inputDescription = QLineEdit(application)
        self.inputDescription.setPlaceholderText("Enter name of the molecule...")
        self.inputDescription.setFixedWidth(250)
        self.submitButton = QPushButton("Add to catalogue", application)
        self.submitButton.clicked.connect(application.onSubmitButtonClicked)
        self.submitButton.setFixedWidth(250)

        self.inputSMILES.setStyleSheet("""
            QLineEdit {
                border: 2px solid #A0A0A0;
                border-radius: 5px;
                padding: 5px;
                margin-bottom: 10px;
            }
        """)

        self.inputDescription.setStyleSheet("""
            QLineEdit {
                border: 2px solid #A0A0A0;
                border-radius: 5px;
                padding: 5px;
                margin-bottom: 10px;
            }
        """)

        self.submitButton.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)

        self.inputVLayout = QVBoxLayout()
        self.inputVLayout.addWidget(self.customMoleculeTitle)
        self.inputVLayout.addWidget(self.inputSMILES)
        self.inputVLayout.addWidget(self.inputDescription)
        self.inputVLayout.addWidget(self.submitButton)

        self.container = QWidget()
        self.container.setLayout(self.inputVLayout)
        self.container.setStyleSheet("""
            padding-top: 0px;
        """)
        # self.container.setFixedWidth(270)
        self.container.setFixedSize(270, 200)

    def getForm(self):
        return self.container

    def getInputSmilesText(self):
        return self.inputSMILES.text()

    def getInputDescriptionText(self):
        return self.inputDescription.text()