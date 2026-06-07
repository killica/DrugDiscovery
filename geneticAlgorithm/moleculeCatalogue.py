from individual import Individual
import json


class MoleculeCatalogue:
    def __init__(self, molecules=None):
        # All available molecules
        self.molecules: list[Individual] = molecules if molecules else []

        # Selected molecules for the first generation
        self.selectedMolecules: list[Individual] = []

        # Newly generated molecules
        self.newGenerationMolecules: list[Individual] = []

    def addMolecule(self, smiles: str, description: str, weights):
        individual = Individual(smiles, description, weights)
        self.molecules.append(individual)
        self.saveToFile(smiles, description)
        return individual

    def saveToFile(self, smiles, description, filepath="../data/molecules.json"):
        try:
            with open(filepath, "r") as file:
                data = json.load(file)
        except FileNotFoundError:
            data = []

        data.append({
            "SMILES": smiles,
            "Description": description
        })

        with open(filepath, "w") as file:
            json.dump(data, file, indent=4)

    def removeFromCatalogue(self, individual: Individual) -> Individual:
        self.molecules.remove(individual)
        return individual

    def addToCatalogue(self, individual: Individual):
        self.molecules.append(individual)

    def removeFromSelected(self, individual: Individual) -> Individual:
        self.selectedMolecules.remove(individual)
        return individual

    def addToSelected(self, individual: Individual):
        self.selectedMolecules.append(individual)

    def selectAll(self):
        self.selectedMolecules.extend(self.molecules)
        self.molecules.clear()

    def sortMolecules(self, weights):
        for ind in (
            self.molecules
            + self.selectedMolecules
            + self.newGenerationMolecules
        ):
            ind.setWeights(weights)

        self.molecules.sort(reverse=True)
        self.selectedMolecules.sort(reverse=True)
        self.newGenerationMolecules.sort(reverse=True)
