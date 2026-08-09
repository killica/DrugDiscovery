from individual import Individual
from catalog_paths import append_to_catalog


class MoleculeCatalogue:
    def __init__(self, molecules=None):
        # All available molecules
        self.molecules: list[Individual] = molecules if molecules else []

        # Selected molecules for the first generation
        self.selectedMolecules: list[Individual] = []

        # Newly generated molecules
        self.newGenerationMolecules: list[Individual] = []

    def addMolecule(self, smiles: str, description: str, weights, fitness_mode=0):
        individual = Individual(smiles, description, weights, fitness_mode=fitness_mode)
        self.molecules.append(individual)
        self.saveToFile(smiles, description)
        return individual

    def saveToFile(self, smiles, description, filepath=None):
        if filepath is not None:
            import json
            try:
                with open(filepath, encoding="utf-8") as file:
                    data = json.load(file)
            except FileNotFoundError:
                data = []
            data.append({"SMILES": smiles, "Description": description})
            with open(filepath, "w", encoding="utf-8") as file:
                json.dump(data, file, indent=4)
                file.write("\n")
            return
        append_to_catalog(smiles, description)

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

    def sortMolecules(self, weights, fitness_mode=None):
        for ind in (
            self.molecules
            + self.selectedMolecules
            + self.newGenerationMolecules
        ):
            if fitness_mode is None:
                ind.setWeights(weights)
            else:
                ind.update_fitness_context(weights, fitness_mode)

        self.molecules.sort(reverse=True)
        self.selectedMolecules.sort(reverse=True)
        self.newGenerationMolecules.sort(reverse=True)
