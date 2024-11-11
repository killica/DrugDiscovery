from fitness import Fitness
from rdkit import Chem

class Individual:
    def __init__(self, smiles, description, weights = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)):
        self.__smiles = smiles
        self.__description = description
        self.__weights = weights
        self.calcFitness()

    def calcFitness(self):
        fit = Fitness(Chem.MolFromSmiles(self.__smiles), self.__weights)
        self.__qed = fit.qed()

    def __lt__(self, other):
        return self.__qed < other.getQED()

    def getSmiles(self):
        return self.__smiles

    def getDescription(self):
        return self.__description

    def getWeights(self):
        return self.__weights

    def getQED(self):
        return self.__qed

    def setWeights(self, weights):
        self.__weights = weights
        self.calcFitness()

    def setDescription(self, description):
        self.__description = description
