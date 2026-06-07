"""Render molecule structure images for reports and analysis."""

from rdkit import Chem
from rdkit.Chem import Draw

try:
    from PyQt5.QtGui import QImage, QPixmap
except ImportError:
    QImage = None
    QPixmap = None

DEFAULT_IMAGE_SIZE = 200


def smiles_to_pixmap(smiles, size=DEFAULT_IMAGE_SIZE):
    if QPixmap is None or not smiles:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol_image = Draw.MolToImage(mol, size=(size, size))
    qimage = QImage(
        mol_image.tobytes(),
        mol_image.width,
        mol_image.height,
        mol_image.width * 3,
        QImage.Format_RGB888,
    )
    return QPixmap.fromImage(qimage)
