from fitness import FITNESS_MODE_QED, Fitness, INVALID_FITNESS
from rdkit import Chem


class Individual:
    def __init__(
        self,
        smiles,
        description,
        weights=(0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95),
        fitness_mode=FITNESS_MODE_QED,
    ):
        self.__smiles = smiles
        self.__description = description
        self.__weights = weights
        self.__fitness_mode = fitness_mode
        self.__fitness = INVALID_FITNESS
        self.calcFitness()

    def calcFitness(self):
        fit = Fitness(
            Chem.MolFromSmiles(self.__smiles),
            self.__weights,
            mode=self.__fitness_mode,
            smiles=self.__smiles,
        )
        self.__fitness = fit.value()

    def update_fitness_context(self, weights=None, fitness_mode=None):
        if weights is not None:
            self.__weights = weights
        if fitness_mode is not None:
            self.__fitness_mode = fitness_mode
        self.calcFitness()

    def __lt__(self, other):
        return self.getFitness() < other.getFitness()

    def getSmiles(self):
        return self.__smiles

    def getDescription(self):
        return self.__description

    def getWeights(self):
        return self.__weights

    def getFitnessMode(self):
        return self.__fitness_mode

    def getFitness(self):
        return self.__fitness

    def getQED(self):
        """Backward-compatible alias; returns the active fitness score."""
        return self.__fitness

    def setFitnessMode(self, fitness_mode):
        if fitness_mode != self.__fitness_mode:
            self.__fitness_mode = fitness_mode
            self.calcFitness()

    def setWeights(self, weights):
        self.__weights = weights
        self.calcFitness()

    def setDescription(self, description):
        self.__description = description

    def setSmiles(self, smiles):
        self.__smiles = smiles
