from PyQt5.QtWidgets import (
    QVBoxLayout,
    QPushButton,
    QLabel,
    QLineEdit,
    QFrame,
    QMessageBox,
)


class NewMoleculeForm:
    """Inline panel for adding a molecule to the catalogue (stage 1 sidebar)."""

    def __init__(self, application):
        self.application = application

        self.outer = QFrame()
        self.outer.setObjectName("addMoleculeSection")
        self.outer.setStyleSheet(
            """
            QFrame#addMoleculeSection {
                border: 1px solid #a5d6a7;
                border-radius: 8px;
                background-color: #f8faf8;
            }
            """
        )

        self.customMoleculeTitle = QLabel("Add molecule to catalogue")
        self.customMoleculeTitle.setStyleSheet(
            "font-size: 15px; font-weight: bold; color: #1b5e20;"
        )
        self.inputSMILES = QLineEdit(self.outer)
        self.inputSMILES.setPlaceholderText("SMILES…")
        self.inputDescription = QLineEdit(self.outer)
        self.inputDescription.setPlaceholderText("Name or description (optional)")
        self.submitButton = QPushButton("Add to catalogue", self.outer)
        self.submitButton.clicked.connect(self._handle_submit)

        self.inputSMILES.setStyleSheet(
            """
            QLineEdit {
                border: 2px solid #a5d6a7;
                border-radius: 8px;
                padding: 8px;
                background-color: #ffffff;
            }
            """
        )
        self.inputDescription.setStyleSheet(
            """
            QLineEdit {
                border: 2px solid #a5d6a7;
                border-radius: 8px;
                padding: 8px;
                background-color: #ffffff;
            }
            """
        )
        self.submitButton.setStyleSheet(
            """
            QPushButton {
                background-color: #2e7d32;
                color: white;
                border: none;
                border-radius: 8px;
                padding: 10px 14px;
                font-size: 14px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #1b5e20; }
            """
        )

        inner = QVBoxLayout(self.outer)
        inner.setContentsMargins(12, 12, 12, 12)
        inner.setSpacing(10)
        inner.addWidget(self.customMoleculeTitle)
        inner.addWidget(self.inputSMILES)
        inner.addWidget(self.inputDescription)
        inner.addWidget(self.submitButton)

    def _handle_submit(self):
        if not self.inputSMILES.text().strip():
            QMessageBox.warning(
                self.outer,
                "Missing SMILES",
                "Enter a SMILES string for the molecule.",
            )
            return
        if self.application.onSubmitButtonClicked():
            self.inputSMILES.clear()
            self.inputDescription.clear()

    def getPanelWidget(self):
        return self.outer

    def getInputSmilesText(self):
        return self.inputSMILES.text()

    def getInputDescriptionText(self):
        return self.inputDescription.text()
