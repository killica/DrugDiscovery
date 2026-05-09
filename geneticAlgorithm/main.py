import sys
import json
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QSlider,
    QDesktopWidget,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QLineEdit,
    QScrollArea,
    QGroupBox,
    QStackedWidget,
    QMessageBox,
    QFrame,
    QSizePolicy,
    QButtonGroup,
    QRadioButton,
)
from PyQt5.QtGui import QImage, QPixmap, QPainter, QColor, QPen, QPalette
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw
from GAConfig import GAConfig
from moleculeBoxes import MoleculeBoxes
from insertMolecule import NewMoleculeForm
from hyperParameters import HyperParameters
from gaParameters import GAParameters
from individual import Individual
from mutationInfo import MutationInfo

STAGE1_WINDOW_SIZE = (1240, 980)
STAGE2_WINDOW_SIZE = (1680, 960)

def apply_light_fusion_theme(app):
    """Use Fusion + a light palette so the UI stays white on macOS/Windows dark mode."""
    app.setStyle("Fusion")
    palette = QPalette()
    window = QColor(255, 255, 255)
    text = QColor(30, 30, 30)
    muted = QColor(120, 120, 120)
    palette.setColor(QPalette.Window, window)
    palette.setColor(QPalette.WindowText, text)
    palette.setColor(QPalette.Base, window)
    palette.setColor(QPalette.AlternateBase, QColor(245, 245, 245))
    palette.setColor(QPalette.ToolTipBase, QColor(255, 255, 220))
    palette.setColor(QPalette.ToolTipText, text)
    palette.setColor(QPalette.Text, text)
    palette.setColor(QPalette.Button, QColor(240, 240, 240))
    palette.setColor(QPalette.ButtonText, text)
    palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
    palette.setColor(QPalette.Link, QColor(0, 100, 200))
    palette.setColor(QPalette.Highlight, QColor(0, 120, 215))
    palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))
    if hasattr(QPalette, "PlaceholderText"):
        palette.setColor(QPalette.PlaceholderText, QColor(127, 127, 127))
    for group in (QPalette.Disabled,):
        palette.setColor(group, QPalette.WindowText, muted)
        palette.setColor(group, QPalette.Text, muted)
        palette.setColor(group, QPalette.ButtonText, muted)
    app.setPalette(palette)


class Application(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Drug Discovery")
        self.resize(800, 600)

        self.molecules = self.readMolecules()
        self.sliderValues = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95]

        self.blockTransfer = False
        self._cancel_evolution = False
        self.fitness_mode_id = 0

        self.moleculeBoxes = MoleculeBoxes(self)
        self.newMoleculeForm = NewMoleculeForm(self)
        self.hyperParamLayout = HyperParameters(self)
        self.gaParameters = GAParameters(self)

        self.mi = MutationInfo()

        self.gaConfig = GAConfig(
            generations=100,
            tournamentSize=4,
            elitismSize=1,
            mutationProbability=0.05,
            rouletteSelection=False,
        )

        self.sbmtBtn = self.newMoleculeForm.submitButton
        self.resBtn = self.hyperParamLayout.resetButton
        self.hyperParamLayout.setSliderTrackWidth(268)

        # --- Stage 1: catalogue + selected (stacked) + sidebar (fitness + add molecule + continue) ---
        self.stage1Sidebar = QWidget()
        self.stage1Sidebar.setMaximumWidth(420)
        self.stage1Sidebar.setMinimumWidth(320)
        self.stage1Sidebar.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        sidebar_layout = QVBoxLayout(self.stage1Sidebar)
        sidebar_layout.setSpacing(12)
        sidebar_layout.setContentsMargins(14, 8, 14, 8)

        fitness_title = QLabel("Fitness")
        fitness_title.setStyleSheet("font-size: 17px; font-weight: bold; color: #1b5e20;")
        sidebar_layout.addWidget(fitness_title)

        self.fitness_mode_group = QButtonGroup(self)
        self.radio_fitness_qed = QRadioButton("QED (weighted)")
        self.radio_fitness_qed.setChecked(True)
        self.radio_fitness_m1 = QRadioButton("Model A (placeholder)")
        self.radio_fitness_m2 = QRadioButton("Model B (placeholder)")
        self.radio_fitness_m3 = QRadioButton("Model C (placeholder)")
        for idx, rb in enumerate(
            (
                self.radio_fitness_qed,
                self.radio_fitness_m1,
                self.radio_fitness_m2,
                self.radio_fitness_m3,
            )
        ):
            self.fitness_mode_group.addButton(rb, idx)
            sidebar_layout.addWidget(rb)
        self.fitness_mode_group.idClicked.connect(self._on_fitness_mode_changed)

        self.fitness_stack = QStackedWidget()
        self.fitness_stack.addWidget(self.hyperParamLayout.getSlidersWidget())
        self.fitness_stack.addWidget(self._fitness_model_placeholder_panel("Model A"))
        self.fitness_stack.addWidget(self._fitness_model_placeholder_panel("Model B"))
        self.fitness_stack.addWidget(self._fitness_model_placeholder_panel("Model C"))
        self.fitness_stack.setMinimumHeight(320)

        self.fitness_stack_scroll = QScrollArea()
        self.fitness_stack_scroll.setWidgetResizable(True)
        self.fitness_stack_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.fitness_stack_scroll.setWidget(self.fitness_stack)
        self.fitness_stack_scroll.setMinimumHeight(400)
        self.fitness_stack_scroll.setMaximumHeight(720)
        self.fitness_stack_scroll.setFrameShape(QFrame.NoFrame)
        self.fitness_stack_scroll.setStyleSheet(
            "QScrollArea { background: transparent; border: 1px solid #c8e6c9; border-radius: 8px; }"
        )
        sidebar_layout.addWidget(self.fitness_stack_scroll, 0)

        sidebar_layout.addSpacing(8)
        sep_actions = QFrame()
        sep_actions.setFrameShape(QFrame.HLine)
        sep_actions.setFixedHeight(1)
        sep_actions.setStyleSheet("background-color: #e0e0e0; border: none;")
        sidebar_layout.addWidget(sep_actions)

        sidebar_layout.addWidget(self.newMoleculeForm.getPanelWidget())
        sidebar_layout.addStretch(1)

        self.continueToNextButton = QPushButton("Continue")
        self.continueToNextButton.setFixedWidth(200)
        self.continueToNextButton.setStyleSheet("""
            QPushButton {
                background-color: #2e7d32;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 12px 16px;
                font-size: 15px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #1b5e20; }
        """)
        self.continueToNextButton.clicked.connect(self.on_continue_to_stage_2)
        sidebar_layout.addWidget(self.continueToNextButton, 0, Qt.AlignHCenter)

        stage1_sep = QFrame()
        stage1_sep.setFrameShape(QFrame.VLine)
        stage1_sep.setFrameShadow(QFrame.Sunken)
        stage1_sep.setStyleSheet("color: #c8c8c8;")
        stage1_sep.setFixedWidth(1)
        stage1_sep.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)

        stage1_row = QHBoxLayout()
        stage1_row.setSpacing(0)
        stage1_row.setContentsMargins(8, 6, 8, 6)
        stage1_row.addWidget(self.moleculeBoxes.getSelectionWidget(), 1)
        stage1_row.addSpacing(20)
        stage1_row.addWidget(stage1_sep, 0, Qt.AlignTop)
        stage1_row.addSpacing(16)
        stage1_row.addWidget(self.stage1Sidebar, 0)

        self.stage1Page = QWidget()
        self.stage1Page.setLayout(stage1_row)

        # --- Stage 2: GA config (left) | evolution views (right); ---
        self.stage2LeftLayout = QVBoxLayout()
        self.stage2LeftLayout.setSpacing(16)
        self.stage2LeftLayout.addWidget(self.gaParameters.getGAParametersWidget())

        self.stage2Left = QWidget()
        self.stage2Left.setLayout(self.stage2LeftLayout)
        self.stage2Left.setMinimumWidth(780)

        self.evolutionColumn = QVBoxLayout()
        self.evolutionColumn.addWidget(self.moleculeBoxes.getPrecedentScrollArea())
        self.evolutionColumn.addWidget(self.moleculeBoxes.getSecondScrollArea())
        self.evolutionColumn.addWidget(self.moleculeBoxes.getBest())

        self.rightWrapper = QWidget()
        self.rightWrapper.setLayout(self.evolutionColumn)
        self.rightWrapper.setFixedHeight(1030)

        stage2_row = QHBoxLayout()
        stage2_row.setSpacing(16)
        stage2_row.addWidget(self.stage2Left, 0, Qt.AlignTop)
        stage2_row.addWidget(self.rightWrapper, 0, Qt.AlignTop)

        self.stage2Page = QWidget()
        self.stage2Page.setLayout(stage2_row)

        self.stack = QStackedWidget()
        self.stack.addWidget(self.stage1Page)
        self.stack.addWidget(self.stage2Page)

        root = QVBoxLayout()
        root.setContentsMargins(8, 8, 8, 8)
        root.addWidget(self.stack)
        self.setLayout(root)
        self.setMinimumSize(1180, 860)
        self.resize(*STAGE1_WINDOW_SIZE)
        self.setAutoFillBackground(True)

        self.stack.setCurrentIndex(0)

        self.show()

    def _fitness_model_placeholder_panel(self, name):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        msg = QLabel(
            f"{name}\n\n"
            "Add controls for this model here (hyperparameters, checkpoint path, "
            "train/eval actions, etc.).\n\n"
            "Ranking still uses QED until this backend is connected."
        )
        msg.setWordWrap(True)
        msg.setStyleSheet("color: #424242; font-size: 13px;")
        lay.addWidget(msg)
        lay.addStretch(1)
        return w

    def _on_fitness_mode_changed(self, mode_id):
        self.fitness_mode_id = mode_id
        self.fitness_stack.setCurrentIndex(mode_id)

    def get_fitness_mode_id(self):
        return self.fitness_mode_id

    def on_continue_to_stage_2(self):
        if len(self.moleculeBoxes.selectedMolecules) == 0:
            QMessageBox.warning(
                self,
                "Nothing selected",
                "Select at least one molecule for the first generation (click cards in the catalogue).",
            )
            return
        self.moleculeBoxes.place_precedent_in_evolution_row()
        self.stack.setCurrentIndex(1)
        self.resize(*STAGE2_WINDOW_SIZE)

    def show_stage_1(self):
        self.moleculeBoxes.place_precedent_in_catalogue_column()
        self.stack.setCurrentIndex(0)
        self.resize(*STAGE1_WINDOW_SIZE)

    def closeEvent(self, event):
        self._cancel_evolution = True
        super().closeEvent(event)

    def paintEvent(self, event):
        super().paintEvent(event)
        if self.stack.currentIndex() != 1:
            return
        painter = QPainter(self)
        painter.setPen(QPen(QColor(128, 128, 128), 2))
        painter.drawLine(10, 685, 800, 685)
        painter.drawLine(800, 20, 800, 880)
        painter.end()

    def onSubmitButtonClicked(self):
        smiles = self.newMoleculeForm.getInputSmilesText().strip()
        description = self.newMoleculeForm.getInputDescriptionText().strip()
        return self.moleculeBoxes.addToCatalogue(smiles, description)

    def readMolecules(self):
        with open("../data/molecules.json", "r") as file:
            data = json.load(file)
        return [Individual(item["SMILES"], item["Description"]) for item in data]

    def getMolecules(self):
        return self.molecules

    def getSliderValues(self):
        return self.sliderValues

    def getMutationInfo(self):
        return self.mi


if __name__ == "__main__":
    app = QApplication(sys.argv)
    apply_light_fusion_theme(app)
    window = Application()
    sys.exit(app.exec_())
