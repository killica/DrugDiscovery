from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QLabel, QLineEdit
from PyQt5.QtCore import Qt

class NewMoleculeForm:
    def __init__(self, application):
        self.customMoleculeTitle = QLabel("Or, add your own molecule:", application)
        self.customMoleculeTitle.setAlignment(Qt.AlignCenter)

        self.customMoleculeTitle.setStyleSheet("font-size: 20px; font-weight: bold; padding-bottom: 10px;")

        self.inputSMILES = QLineEdit(application)
        self.inputSMILES.setPlaceholderText("Enter SMILES format of designated molecule...")
        self.inputDescription = QLineEdit(application)
        self.inputDescription.setPlaceholderText("Enter name of the molecule...")
        self.submitButton = QPushButton("Add to catalogue", application)
        self.submitButton.clicked.connect(application.onSubmitButtonClicked)

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

    def getLayout(self):
        return self.inputVLayout

    def getSubmitButton(self):
        return self.submitButton

    def getInputSmilesText(self):
        return self.inputSMILES.text()

    def getInputDescriptionText(self):
        return self.inputDescription.text()