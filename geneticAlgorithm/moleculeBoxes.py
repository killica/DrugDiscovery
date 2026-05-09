import json
import copy
import re
import random
from PyQt5.QtWidgets import (
    QGroupBox,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QGraphicsDropShadowEffect,
    QWidget,
    QGridLayout,
    QScrollArea,
    QPushButton,
    QProgressBar,
    QApplication,
    QFrame,
    QSizePolicy,
    QDialog,
    QDialogButtonBox,
    QLineEdit,
)
from PyQt5.QtGui import QImage, QPixmap, QColor, QFont, QPainter, QPainterPath, QPen, QCursor
from PyQt5.QtCore import Qt, QRectF, QEvent, QTimer

MOLECULE_IMAGE_SIZE = 200
MOLECULE_IMAGE_CORNER_RADIUS = 14

STAGE1_SCROLL_WIDTH = 760
STAGE1_SCROLL_HEIGHT = 420

# Distinct framed scroll panes + separator so the two stage-1 lists are obviously separate.
STAGE1_SCROLL_QSS_CATALOGUE = """
    QScrollArea {
        border: 1px solid #b0bec5;
        border-radius: 10px;
        background-color: #fafbfc;
    }
"""
STAGE1_SCROLL_QSS_SELECTED = """
    QScrollArea {
        border: 1px solid #43a047;
        border-radius: 10px;
        background-color: #f3faf4;
    }
"""

MOLECULE_BOX_QSS_NORMAL = """
    QGroupBox {
        background-color: rgb(255, 255, 255);
        border: 2px solid green;
        border-radius: 15px;
        padding: 10px;
    }
"""
MOLECULE_BOX_QSS_HOVER = """
    QGroupBox {
        background-color: #e8f5e9;
        border: 2px solid green;
        border-radius: 15px;
        padding: 10px;
    }
"""

class RoundedMoleculeImage(QLabel):
    """Renders a square pixmap clipped to rounded corners."""

    def __init__(self, pixmap, parent=None):
        super().__init__(parent)
        self._pixmap = pixmap
        self.setFixedSize(MOLECULE_IMAGE_SIZE, MOLECULE_IMAGE_SIZE)
        self.setAttribute(Qt.WA_TranslucentBackground, True)

    def setMoleculePixmap(self, pixmap):
        self._pixmap = pixmap
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        w, h = float(self.width()), float(self.height())
        r = float(MOLECULE_IMAGE_CORNER_RADIUS)
        path = QPainterPath()
        path.addRoundedRect(QRectF(0, 0, w, h), r, r)
        painter.setClipPath(path)

        if self._pixmap is None or self._pixmap.isNull():
            painter.fillPath(path, QColor(245, 245, 245))
        else:
            scaled = self._pixmap.scaled(
                self.size(),
                Qt.IgnoreAspectRatio,
                Qt.SmoothTransformation,
            )
            painter.drawPixmap(0, 0, scaled)

        painter.setClipping(False)
        painter.setPen(QPen(QColor(218, 222, 226), 1))
        painter.setBrush(Qt.NoBrush)
        painter.drawRoundedRect(QRectF(0.5, 0.5, w - 1, h - 1), r, r)
        painter.end()
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors
from rdkit.DataStructs import FingerprintSimilarity
from fitness import Fitness
from individual import Individual
import geneticAlgorithm
from mutationInfo import MutationInfo
from datetime import datetime

class ClickableGroupBox(QGroupBox):
    def __init__(self, moleculeBoxes, index, ind, parent=None):
        super().__init__(parent)
        self.moleculeBoxes = moleculeBoxes
        self.index = index
        self.ind = ind
        self.description_label = None
        self._molecule_hover_active = False
        self.setAttribute(Qt.WA_StyledBackground, True)

    def _install_descendant_hover_tracking(self):
        self.installEventFilter(self)
        for w in self.findChildren(QWidget):
            w.installEventFilter(self)

    def eventFilter(self, watched, event):
        et = event.type()
        if et == QEvent.Enter:
            self._set_molecule_hover(True)
        elif et == QEvent.Leave:
            QTimer.singleShot(0, self._sync_molecule_hover_from_cursor)
        return False

    def _sync_molecule_hover_from_cursor(self):
        w = QApplication.widgetAt(QCursor.pos())
        if w is not None and (w is self or self.isAncestorOf(w)):
            self._set_molecule_hover(True)
        else:
            self._set_molecule_hover(False)

    def _set_molecule_hover(self, active):
        if self._molecule_hover_active == active:
            return
        self._molecule_hover_active = active
        self.setStyleSheet(MOLECULE_BOX_QSS_HOVER if active else MOLECULE_BOX_QSS_NORMAL)
        if self.description_label is None:
            return
        if active and not self.moleculeBoxes.application.blockTransfer:
            self.description_label.setStyleSheet(
                "color: darkgreen; font-weight: bold; background: transparent;"
            )
        else:
            self.description_label.setStyleSheet(
                "color: #111111; font-weight: normal; background: transparent;"
            )

    def mousePressEvent(self, event):
        if self.moleculeBoxes.application.blockTransfer:
            return
        app = self.moleculeBoxes.application
        if hasattr(app, "stack") and app.stack.currentIndex() == 2:
            super().mousePressEvent(event)
            return
        if self.ind != 0 and self.ind != 1:
            return
        if event.button() == Qt.LeftButton:
            self.moleculeBoxes.removeBoxes()
            self.moleculeBoxes.removeSelectedBoxes()
            if self.ind == 0:
                self.moleculeBoxes.selectedMolecules.append(self.moleculeBoxes.molecules[self.index])
                self.moleculeBoxes.molecules.pop(self.index)
            elif self.ind == 1:
                self.moleculeBoxes.molecules.append(self.moleculeBoxes.selectedMolecules[self.index])
                self.moleculeBoxes.selectedMolecules.pop(self.index)
            self.moleculeBoxes.loadBoxes(tuple(self.moleculeBoxes.application.getSliderValues()))
            self.moleculeBoxes.loadSelectedBoxes(tuple(self.moleculeBoxes.application.getSliderValues()))

        super().mousePressEvent(event)

class MoleculeBoxes(QWidget):
    def __init__(self, application):
        super().__init__(application)
        self.application = application
        self.molecules = application.molecules
        self.windowWidth = application.width()
        self.boxWidth = 220
        self.columnsPerRow = 3
        # Avoid re-running RDKit MolToImage for the same SMILES on every layout rebuild (slider / transfer).
        self._structure_pixmap_cache = {}
        self.layout = QGridLayout()
        self.layout.setContentsMargins(6, 6, 6, 6)

        self.selectedMolecules = []
        self.newGenerationMolecules = []
        self.newGenerationBoxes = []

        self.selectAllButton = QPushButton("Select all")
        self.selectAllButton.setFixedWidth(100)
        self.selectAllButton.setStyleSheet("""
            QPushButton {
                background-color: #696969;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
            }
            QPushButton:hover {
                background-color: #404040;
            }
        """)
        self.selectAllButton.clicked.connect(self.onSelectAllButtonClicked)

        self.sampleButton = QPushButton("Sample")
        self.sampleButton.setFixedWidth(88)
        self.sampleButton.setStyleSheet("""
            QPushButton {
                background-color: #5c8a8a;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
            }
            QPushButton:hover {
                background-color: #3d6666;
            }
        """)
        self.sampleButton.clicked.connect(self.onSampleButtonClicked)

        self.removeAllButton = QPushButton("Remove all")
        self.removeAllButton.setFixedWidth(108)
        self.removeAllButton.setStyleSheet("""
            QPushButton {
                background-color: #8d4d4d;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
            }
            QPushButton:hover {
                background-color: #6d3333;
            }
        """)
        self.removeAllButton.clicked.connect(self.onRemoveAllButtonClicked)

        self.catalogueLabel = QLabel("Catalogue")
        self.catalogueLabel.setStyleSheet(
            "font-size: 18px; font-weight: bold; margin: 0; padding: 0 0 2px 0;"
        )

        self.scrollArea = QScrollArea()
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setFixedSize(STAGE1_SCROLL_WIDTH, STAGE1_SCROLL_HEIGHT)
        self.scrollArea.setStyleSheet(STAGE1_SCROLL_QSS_CATALOGUE)

        self.catalogueGridHost = QWidget()
        self.catalogueGridHost.setLayout(self.layout)
        self.scrollWidget = QWidget()
        self.catalogueOuterLayout = QVBoxLayout(self.scrollWidget)
        self.catalogueOuterLayout.setContentsMargins(0, 0, 0, 0)
        self.catalogueOuterLayout.setSpacing(0)
        self.catalogueOuterLayout.addWidget(self.catalogueGridHost, 0, Qt.AlignTop)
        self.catalogueOuterLayout.addStretch(1)
        self.scrollArea.setWidget(self.scrollWidget)

        self.selectedSectionLabel = QLabel("Selected for first generation")
        self.selectedSectionLabel.setStyleSheet(
            "font-size: 18px; font-weight: bold; color: darkgreen; margin: 0; padding: 2px 0 0 0;"
        )

        self.precedentLayout = QGridLayout()
        self.precedentLayout.setContentsMargins(6, 6, 6, 6)

        self.precedentScrollArea = QScrollArea()
        self.precedentScrollArea.setWidgetResizable(True)
        self.precedentScrollArea.setFixedSize(STAGE1_SCROLL_WIDTH, STAGE1_SCROLL_HEIGHT)
        self.precedentScrollArea.setStyleSheet(STAGE1_SCROLL_QSS_SELECTED)

        self.precedentGridHost = QWidget()
        self.precedentGridHost.setLayout(self.precedentLayout)
        self.precedentScrollWidget = QWidget()
        self.precedentOuterLayout = QVBoxLayout(self.precedentScrollWidget)
        self.precedentOuterLayout.setContentsMargins(0, 0, 0, 0)
        self.precedentOuterLayout.setSpacing(0)
        self.precedentOuterLayout.addWidget(self.precedentGridHost, 0, Qt.AlignTop)
        self.precedentOuterLayout.addStretch(1)
        self.precedentScrollWidget.setMinimumWidth(0)
        self.precedentScrollArea.setWidget(self.precedentScrollWidget)

        self.precedentLabel = QLabel("1. generation")
        self.precedentLabel.setStyleSheet("font-size: 20px; font-weight: bold; color: darkgreen;")

        self.stage1SectionSeparator = QFrame()
        self.stage1SectionSeparator.setObjectName("stage1SectionSeparator")
        self.stage1SectionSeparator.setFrameShape(QFrame.NoFrame)
        self.stage1SectionSeparator.setFixedHeight(2)
        self.stage1SectionSeparator.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.stage1SectionSeparator.setStyleSheet(
            "QFrame#stage1SectionSeparator { background-color: #81c784; border-radius: 1px; }"
        )

        self.catalogueHeaderRow = QWidget()
        _cat_hdr = QHBoxLayout(self.catalogueHeaderRow)
        _cat_hdr.setContentsMargins(0, 0, 0, 0)
        _cat_hdr.setSpacing(10)
        _cat_hdr.addWidget(self.catalogueLabel, 0, Qt.AlignVCenter)
        _cat_hdr.addStretch(1)
        _cat_hdr.addWidget(self.selectAllButton, 0, Qt.AlignVCenter)
        _cat_hdr.addWidget(self.sampleButton, 0, Qt.AlignVCenter)

        self.selectedHeaderRow = QWidget()
        _sel_hdr = QHBoxLayout(self.selectedHeaderRow)
        _sel_hdr.setContentsMargins(0, 0, 0, 0)
        _sel_hdr.setSpacing(10)
        _sel_hdr.addWidget(self.selectedSectionLabel, 0, Qt.AlignVCenter)
        _sel_hdr.addStretch(1)
        _sel_hdr.addWidget(self.removeAllButton, 0, Qt.AlignVCenter)

        self.catalogueVBox = QVBoxLayout()
        self.catalogueVBox.setSpacing(4)
        self.catalogueVBox.setContentsMargins(0, 0, 4, 0)
        self.catalogueVBox.addWidget(self.catalogueHeaderRow)
        self.catalogueVBox.addWidget(self.scrollArea)
        self.catalogueVBox.addSpacing(12)
        self.catalogueVBox.addWidget(self.stage1SectionSeparator)
        self.catalogueVBox.addSpacing(10)
        self.catalogueVBox.addWidget(self.selectedHeaderRow)
        self.catalogueVBox.addWidget(self.precedentScrollArea)

        self.container = QWidget()
        self.container.setLayout(self.catalogueVBox)
        self.container.setMinimumWidth(STAGE1_SCROLL_WIDTH)
        self.container.setMinimumSize(STAGE1_SCROLL_WIDTH, STAGE1_SCROLL_HEIGHT * 2 + 112)
        self.precedentLabel.hide()

        self.evolutionGen1VBox = QVBoxLayout()
        self.evolutionGen1VBox.setSpacing(8)
        self.evolutionGen1VBox.setContentsMargins(0, 0, 0, 0)
        self.rightCont1 = QWidget()
        self.rightCont1.setLayout(self.evolutionGen1VBox)
        self.rightCont1.setMinimumSize(200, 160)
        self.rightCont1.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.secondLayout = QGridLayout()
        self.secondLayout.setContentsMargins(0, 0, 0, 0)
        self.secondLayout.setHorizontalSpacing(12)
        self.secondLayout.setVerticalSpacing(12)

        self.secondScrollArea = QScrollArea()
        self.secondScrollArea.setWidgetResizable(True)
        self.secondScrollArea.setMinimumSize(0, 160)
        self.secondScrollArea.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.secondScrollWidget = QWidget()
        self.secondScrollWidget.setLayout(self.secondLayout)
        self.secondScrollWidget.setMinimumWidth(0)
        self.secondScrollArea.setWidget(self.secondScrollWidget)

        self.secondLabel = QLabel("2. generation")
        self.secondLabel.setStyleSheet("font-size: 20px; font-weight: bold; color: darkgreen;")

        self.evolutionGen2VBox = QVBoxLayout()
        self.evolutionGen2VBox.setSpacing(8)
        self.evolutionGen2VBox.setContentsMargins(0, 0, 0, 0)
        self.evolutionGen2VBox.addWidget(self.secondLabel, 0, Qt.AlignLeft)
        self.evolutionGen2VBox.addWidget(self.secondScrollArea, 1)

        self.rightCont2 = QWidget()
        self.rightCont2.setLayout(self.evolutionGen2VBox)
        self.rightCont2.setMinimumSize(200, 160)
        self.rightCont2.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.evolutionControlsLayout = QVBoxLayout()
        self.evolutionControlsLayout.setSpacing(0)
        self.evolutionControlsLayout.setContentsMargins(0, 4, 0, 4)
        self.rightCont3 = QWidget()

        self.loadBoxes()
        self.loadSelectedBoxes()
       
    def loadBoxes(self, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self.boxes = []
        for individual in self.molecules:
            individual.setWeights(weights)
        self.molecules.sort(reverse=True)
        for index, individual in enumerate(self.molecules):
            row = index // 3
            col = index % 3
            smiles = individual.getSmiles()
            description = individual.getDescription()
            qed = individual.getQED()
            self.moleculeBox = self.createMoleculeBox(smiles, description, qed, index, 0)
            self.boxes.append(self.moleculeBox)
            self.layout.addWidget(self.moleculeBox, row, col)

    def loadSelectedBoxes(self, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self.selectedBoxes = []
        for individual in self.selectedMolecules:
            individual.setWeights(weights)
        self.selectedMolecules.sort(reverse=True)
        for index, individual in enumerate(self.selectedMolecules):
            row = index // 3
            col = index % 3
            smiles = individual.getSmiles()
            description = individual.getDescription()
            qed = individual.getQED()
            self.selectedMoleculeBox = self.createMoleculeBox(smiles, description, qed, index, 1)
            self.selectedBoxes.append(self.selectedMoleculeBox)
            self.precedentLayout.addWidget(self.selectedMoleculeBox, row, col)

    def _clear_second_generation_grid(self):
        """Remove every widget from the 2nd-generation layout without relying on
        newGenerationBoxes."""
        while self.secondLayout.count():
            item = self.secondLayout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        self.newGenerationBoxes = []

    def loadNewGeneration(self, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self._clear_second_generation_grid()
        QApplication.processEvents()
        if len(self.newGenerationMolecules) > 0:
            self.newGenerationMolecules.sort(reverse=True)
            mols = self.newGenerationMolecules
            n = len(mols)
            cpr = self.columnsPerRow
            for index in range(n):
                row = index // cpr
                col = index % cpr
                individual = mols[index]
                smiles = individual.getSmiles()
                description = individual.getDescription()
                individual.setWeights(weights)
                qed = individual.getQED()
                box = self.createMoleculeBox(smiles, description, qed, index, -1)
                self.newGenerationBoxes.append(box)
                self.secondLayout.addWidget(box, row, col)

            best_idx = -1
            if getattr(self, "bestBox", None) is not None:
                best_idx = self.evolutionControlsLayout.indexOf(self.bestBox)
                self.evolutionControlsLayout.removeWidget(self.bestBox)
                self.bestBox.deleteLater()
            self.bestBox = self.createMoleculeBox(
                self.newGenerationMolecules[0].getSmiles(),
                "Current best",
                self.newGenerationMolecules[0].getQED(),
                0,
                -1,
            )
            if best_idx >= 0:
                self.evolutionControlsLayout.insertWidget(best_idx, self.bestBox, 0, Qt.AlignTop)
            else:
                pc_idx = self.evolutionControlsLayout.indexOf(getattr(self, "progressCnt", None))
                if pc_idx >= 0:
                    self.evolutionControlsLayout.insertWidget(pc_idx, self.bestBox, 0, Qt.AlignTop)
                else:
                    self.evolutionControlsLayout.addWidget(self.bestBox, 0, Qt.AlignTop)
            self.secondScrollWidget.updateGeometry()
            QApplication.processEvents()
       
    def removeBoxes(self):
        for box in self.boxes:
            self.layout.removeWidget(box)
            box.deleteLater()
        self.boxes = []

    def removeSelectedBoxes(self):
        for selectedBox in self.selectedBoxes:
            self.precedentLayout.removeWidget(selectedBox)
            selectedBox.deleteLater()
        self.selectedBoxes = []

    def removeNewGenerationBoxes(self):
        self._clear_second_generation_grid()

    def getSelectionWidget(self):
        return self.container

    def getPrecedentScrollArea(self):
        return self.rightCont1

    def place_precedent_in_evolution_row(self):
        """Move the 1st-generation scroll from the catalogue column to the evolution header row (stage 2)."""
        if self.precedentScrollArea.parentWidget() == self.rightCont1:
            return
        self.catalogueVBox.removeWidget(self.precedentScrollArea)
        self.precedentScrollArea.setParent(None)
        while self.evolutionGen1VBox.count():
            item = self.evolutionGen1VBox.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
        self.evolutionGen1VBox.addWidget(self.precedentLabel, 0, Qt.AlignLeft)
        self.evolutionGen1VBox.addWidget(self.precedentScrollArea, 1)
        self.precedentLabel.show()
        self.precedentScrollArea.show()
        self.precedentScrollArea.setMinimumSize(0, 160)
        self.precedentScrollArea.setMaximumSize(16777215, 16777215)
        self.precedentScrollArea.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def place_precedent_in_catalogue_column(self):
        """Put the selected-molecules scroll back under the catalogue (stage 1)."""
        if self.precedentScrollArea.parentWidget() == self.container:
            return
        while self.evolutionGen1VBox.count():
            item = self.evolutionGen1VBox.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
        self.catalogueVBox.addWidget(self.precedentScrollArea)
        self.precedentLabel.hide()
        self.precedentScrollArea.show()
        self.precedentScrollArea.setFixedSize(STAGE1_SCROLL_WIDTH, STAGE1_SCROLL_HEIGHT)
        self.precedentScrollArea.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

    def getSecondScrollArea(self):
        return self.rightCont2

    def getBest(self):
        return self.rightCont3

    def _structure_pixmap_for_smiles(self, smiles):
        if not smiles:
            return None
        if smiles in self._structure_pixmap_cache:
            return self._structure_pixmap_cache[smiles]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        molImage = Draw.MolToImage(mol, size=(MOLECULE_IMAGE_SIZE, MOLECULE_IMAGE_SIZE))
        qimage = QImage(molImage.tobytes(), molImage.width, molImage.height, molImage.width * 3, QImage.Format_RGB888)
        pm = QPixmap.fromImage(qimage)
        self._structure_pixmap_cache[smiles] = pm
        return pm

    def refreshQEDForAll(self, weights):

        self.removeBoxes()
        self.removeSelectedBoxes()
        self.loadBoxes(weights)
        self.loadSelectedBoxes(weights)

    def createMoleculeBox(self, smiles, description, qed, index, ind):
        box = ClickableGroupBox(self, index, ind)
        box.setFixedWidth(230)
        boxLayout = QVBoxLayout()
        boxLayout.setSpacing(6)

        pm = self._structure_pixmap_for_smiles(smiles)
        if pm is not None:
            imageLabel = RoundedMoleculeImage(pm)
        else:
            imageLabel = RoundedMoleculeImage(None)

        # Shadow on the image frame only (not the whole QGroupBox — keeps text labels reliable).
        imageFrame = QWidget()
        imageFrame.setFixedSize(MOLECULE_IMAGE_SIZE, MOLECULE_IMAGE_SIZE)
        imageFrameLayout = QVBoxLayout(imageFrame)
        imageFrameLayout.setContentsMargins(0, 0, 0, 0)
        imageFrameLayout.addWidget(imageLabel)
        shadowEffect = QGraphicsDropShadowEffect()
        shadowEffect.setOffset(3, 5)
        shadowEffect.setBlurRadius(20)
        shadowEffect.setColor(QColor(0, 0, 0, 95))
        imageFrame.setGraphicsEffect(shadowEffect)

        imageRow = QHBoxLayout()
        imageRow.setContentsMargins(0, 0, 0, 0)
        imageRow.addStretch(1)
        imageRow.addWidget(imageFrame, 0, Qt.AlignHCenter)
        imageRow.addStretch(1)

        descriptionLabel = QLabel(description + "\n" + "QED: " + str(round(qed, 4)))
        descriptionLabel.setWordWrap(True)
        descriptionLabel.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        descriptionLabel.setStyleSheet(
            "color: #111111; font-weight: normal; background: transparent;"
        )
        box.description_label = descriptionLabel

        boxLayout.setContentsMargins(0, 0, 0, 0)
        boxLayout.addLayout(imageRow)
        boxLayout.addWidget(descriptionLabel)
        box.setLayout(boxLayout)

        box.setStyleSheet(MOLECULE_BOX_QSS_NORMAL)
        box._install_descendant_hover_tracking()
        return box
   
    def onSelectAllButtonClicked(self):
        if self.application.blockTransfer:
            return
        self.removeBoxes()
        self.removeSelectedBoxes()
        for i in range(len(self.molecules)):
            self.selectedMolecules.append(self.molecules[i])
            # self.molecules.pop(i)
        self.molecules = []
        self.loadBoxes(tuple(self.application.getSliderValues()))
        self.loadSelectedBoxes(tuple(self.application.getSliderValues()))

    def onSampleButtonClicked(self):
        if self.application.blockTransfer:
            return
        n_cat = len(self.molecules)
        if n_cat == 0:
            QMessageBox.information(
                self,
                "Sample",
                "There are no molecules in the catalogue to sample.",
            )
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Random sample")
        dlg.setModal(True)
        layout = QVBoxLayout(dlg)
        info = QLabel(
            f"The catalogue has {n_cat} molecule(s). "
            "Choose how many to select at random (each picked at most once)."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        initial = min(5, n_cat)
        dlg._last_good = initial
        dlg._chosen_k = initial

        count_edit = QLineEdit(dlg)
        count_edit.setText(str(initial))
        count_edit.setPlaceholderText(f"1–{n_cat}")

        def _sync_last_good_from_field():
            t = count_edit.text().strip()
            if not t:
                return
            try:
                v = int(t)
            except ValueError:
                return
            if 1 <= v <= n_cat:
                dlg._last_good = v

        count_edit.editingFinished.connect(_sync_last_good_from_field)

        def try_ok():
            t = count_edit.text().strip()
            if not t:
                QMessageBox.warning(
                    dlg,
                    "Invalid number",
                    f"Enter a whole number between 1 and {n_cat}.",
                )
                count_edit.setText(str(dlg._last_good))
                return
            try:
                v = int(t)
            except ValueError:
                QMessageBox.warning(
                    dlg,
                    "Invalid number",
                    f"'{t}' is not a valid whole number. Use an integer from 1 to {n_cat}.",
                )
                count_edit.setText(str(dlg._last_good))
                return
            if v < 1 or v > n_cat:
                QMessageBox.warning(
                    dlg,
                    "Invalid number",
                    f"Enter a number between 1 and {n_cat} (you entered {v}).",
                )
                count_edit.setText(str(dlg._last_good))
                return
            dlg._chosen_k = v
            dlg.accept()

        layout.addWidget(count_edit)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(try_ok)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        if dlg.exec_() != QDialog.Accepted:
            return

        k = min(dlg._chosen_k, len(self.molecules))
        if not self.molecules or k < 1:
            return

        self.removeBoxes()
        self.removeSelectedBoxes()

        n = len(self.molecules)
        pick_idx = set(random.sample(range(n), k))
        to_move = [self.molecules[i] for i in range(n) if i in pick_idx]
        self.molecules = [self.molecules[i] for i in range(n) if i not in pick_idx]
        self.selectedMolecules.extend(to_move)

        self.loadBoxes(tuple(self.application.getSliderValues()))
        self.loadSelectedBoxes(tuple(self.application.getSliderValues()))

    def onRemoveAllButtonClicked(self):
        if self.application.blockTransfer:
            return
        if not self.selectedMolecules:
            return
        self.removeBoxes()
        self.removeSelectedBoxes()
        self.molecules.extend(self.selectedMolecules)
        self.selectedMolecules = []
        self.loadBoxes(tuple(self.application.getSliderValues()))
        self.loadSelectedBoxes(tuple(self.application.getSliderValues()))

    def addToCatalogue(self, smiles, description):
        molecule = None
        try:
            molecule = Chem.MolFromSmiles(smiles)
        except ValueError:
            return False
        if molecule is None:
            return False
        if not description:
            description = "Unknown"

        self.removeBoxes()
        self.removeSelectedBoxes()
        self.molecules.append(Individual(smiles, description, self.application.getSliderValues()))
        with open('../data/molecules.json', 'r') as file:
            data = json.load(file)
        data.append({
            "SMILES": smiles,
            "Description": description
        })
        with open('../data/molecules.json', 'w') as file:
            json.dump(data, file, indent=4)
        self.loadBoxes()
        self.loadSelectedBoxes()
        return True

    def onGenerateButtonClicked(self):
        self.saveLabel.setStyleSheet("color: transparent; font-style: italic;")
        self.removeSelectedBoxes()
        self.selectedMolecules = []
        for ind in self.newGenerationMolecules:
            self.selectedMolecules.append(Individual(ind.getSmiles(), ind.getDescription(), tuple(self.application.getSliderValues())))
        self.loadSelectedBoxes(tuple(self.application.getSliderValues()))

        n = len(self.selectedMolecules)
        self.individualProgress.setMaximum(max(n, 1))
        self.individualProgress.setValue(0)
        self.individualLabel.setText(f"Individual: 0/{n}")
        QApplication.processEvents()

        self.newGenerationMolecules = geneticAlgorithm.geneticAlgorithm(
            self.selectedMolecules,
            True,
            self.application.gaConfig.generations,
            self.application.gaConfig.rouletteSelection,
            self.application.gaConfig.tournamentSize,
            self.application.gaConfig.elitismSize,
            self.application.gaConfig.mutationProbability,
            self.application.getMutationInfo(),
            self.individualLabel,
            self.individualProgress,
            cancel_check=lambda: getattr(self.application, "_cancel_evolution", False),
        )

        if getattr(self.application, "_cancel_evolution", False):
            return

        self.loadNewGeneration(tuple(self.application.getSliderValues()))
        labelText = self.secondLabel.text()
        self.precedentLabel.setText(labelText)
        # Regular expression to match a number at the start of the string
        match = re.match(r'^\d+', labelText)
        # Check if a match was found and extract the number
        labelNumber = int(match.group(0)) + 1
        self.secondLabel.setText(str(labelNumber) + ". generation")
        self.generationLabel.setText(f"Generation: {labelNumber}/{self.application.gaConfig.generations}")
        self.generationProgress.setValue(labelNumber)

    def onFinalButtonClicked(self):
        self.application._cancel_evolution = False
        self.saveLabel.setStyleSheet("color: transparent; font-style: italic;")
        labelText = self.secondLabel.text()
        # Regular expression to match a number at the start of the string
        match = re.match(r'^\d+', labelText)
        # Check if a match was found and extract the number
        labelNumber = int(match.group(0))
        for _ in range(self.application.gaConfig.generations - labelNumber):
            if getattr(self.application, "_cancel_evolution", False):
                break
            self.onGenerateButtonClicked()

    def onSaveButtonClicked(self):
        formattedDatetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open('results/best_candidate_molecules.txt', 'a') as candidatesFile:
            candidatesFile.write(f"SMILES: {self.newGenerationMolecules[0].getSmiles()}\nQED: {round(self.newGenerationMolecules[0].getQED(), 4)}\nDate created: {formattedDatetime}\nParameter weights:")
            for value in list(self.newGenerationMolecules[0].getWeights()):
                candidatesFile.write(f"{value} ")
            candidatesFile.write("\n-------------------------------------------\n")
        self.saveLabel.setStyleSheet("color: green; font-style: italic;")

        # Create a QMessageBox
        msgBox = QMessageBox(self)
        # Set the icon for the dialog
        msgBox.setIcon(QMessageBox.Question)
        # Set the window title
        msgBox.setWindowTitle("Population diversity")
        # Set the message in the dialog
        msgBox.setText("Do you want to calculate Tanimoto similarity coefficient for current generation?")
        # Add Yes and No buttons
        msgBox.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        # Show the message box and capture the response
        response = msgBox.exec_()
        if response == QMessageBox.Yes:
            self.tanimoto() 
        
        # Create a QMessageBox
        msgBoxAvg = QMessageBox(self)
        # Set the icon for the dialog
        msgBoxAvg.setIcon(QMessageBox.Question)
        # Set the window title
        msgBoxAvg.setWindowTitle("Average QED coefficient")
        # Set the message in the dialog
        msgBoxAvg.setText("Do you want to average QED coefficient for current generation?")
        # Add Yes and No buttons
        msgBoxAvg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        # Show the message box and capture the response
        responseAvg = msgBoxAvg.exec_()
        if responseAvg == QMessageBox.Yes:
            self.calculateAverageQED()

    def tanimoto(self):
        smilesList = [s.getSmiles() for s in self.newGenerationMolecules]
        mols = [Chem.MolFromSmiles(s) for s in smilesList]
        fingerprints = [rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits = 2048) for mol in mols]
        # Calculate pairwise Tanimoto similarity
        similarities = []
        for i in range(len(fingerprints)):
            for j in range(i+1, len(fingerprints)):
                similarity = FingerprintSimilarity(fingerprints[i], fingerprints[j])
                similarities.append(similarity)
        with open('results/tanimoto.txt', 'a') as tanimotoFile:
            tanimotoFile.write("[")
            for sim in similarities:
                tanimotoFile.write(f"{sim}, ")
            tanimotoFile.write("]\n------------------------------------\n")

    def calculateAverageQED(self):
        formattedDatetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        coeffList = [s.getQED() for s in self.newGenerationMolecules]
        avgQED = sum(coeffList) / len(coeffList)
        with open('results/averageQED.txt', 'a') as tanimotoFile:
            tanimotoFile.write(f"{avgQED}\n{formattedDatetime}\n------------------------------------\n")

    def onRestartButtonClicked(self):
        self.saveLabel.setStyleSheet("color: transparent; font-style: italic;")
        self.generationLabel.setText(f"Generation: 1/{self.application.gaConfig.generations}")
        self.generationProgress.setValue(1)
        self.individualLabel.setText(f"Individual: 0/{len(self.selectedMolecules)}")
        self.individualProgress.setValue(0)
        self.removeBoxes()
        self.removeSelectedBoxes()
        self.removeNewGenerationBoxes()
        self.molecules = self.application.readMolecules()
        self.loadBoxes(self.application.getSliderValues())
        self.application.molecules = []
        self.selectedMolecules = []
        self.newGenerationMolecules = []
        self.bestBox.deleteLater()
        self.generateButton.setDisabled(True)
        self.generateButton.setStyleSheet("Color: #757575;")
        self.finalButton.setDisabled(True)
        self.finalButton.setStyleSheet("Color: #757575;")
        self.saveButton.setDisabled(True)
        self.saveButton.setStyleSheet("Color: #757575;")
        self.restartButton.setDisabled(True)
        self.restartButton.setStyleSheet("Color: #757575;")
        self.loadNewGeneration()
        self.precedentLabel.setText("1. generation")
        self.secondLabel.setText("2. generation")
        self.application.gaParameters.launchButton.setDisabled(False)
        self.application.gaParameters.launchButton.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 16px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)
        self.application.gaParameters.rouletteCheckBox.setDisabled(False)
        self.application.gaParameters.rouletteCheckBox.setStyleSheet("""
            QCheckBox {
                text-decoration: none;
            }
            QCheckBox::indicator {
                width: 20px;
                height: 20px;
                border: 2px solid #777;
                border-radius: 5px;
                background-color: white;
            }
            QCheckBox::indicator:checked {
                background-color: lightgray;
                border: 2px solid gray;
            }
            QCheckBox::indicator:unchecked {
                background-color: white;
                border: 2px solid gray;
            }
        """)
        self.application.gaParameters.generationSpin.setDisabled(False)
        self.application.gaParameters.tournamentSpin.setDisabled(False)
        self.application.gaParameters.elitismSpin.setDisabled(False)
        self.application.gaParameters.mutationLineEdit.setDisabled(False)
        self.application.sbmtBtn.setDisabled(False)
        self.application.resBtn.setDisabled(False)
        self.application.blockTransfer = False
        while self.evolutionControlsLayout.count():
            item = self.evolutionControlsLayout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
        self.application.show_stage_1()