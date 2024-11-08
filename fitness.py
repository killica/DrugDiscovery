from rdkit import Chem
from rdkit.Chem import QED

class Fitness:
    def __init__(self, molecule, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self.molecule = molecule
        self.weights = weights

    def qed(self):
        return QED.qed(self.molecule, self.weights)