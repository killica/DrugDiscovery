from enum import Enum


class CrossoverMode(Enum):
    SMILES = 0
    SELFIES = 1
    GRAPH = 2
    BRICS = 3


class MutationMode(Enum):
    SMILES = 0
    SELFIES = 1
    GRAPH = 2
    BRICS = 3


class GAConfig:
    def __init__(
        self,
        generations: int = 100,
        tournamentSize: int = 4,
        elitismSize: int = 1,
        mutationProbability: float = 0.05,
        rouletteSelection: bool = False,
        crossoverMode: CrossoverMode = CrossoverMode.SELFIES,
        mutationMode: MutationMode = MutationMode.SMILES,
    ):
        self.generations = generations
        self.tournamentSize = tournamentSize
        self.elitismSize = elitismSize
        self.mutationProbability = mutationProbability
        self.rouletteSelection = rouletteSelection
        self.crossoverMode = crossoverMode
        self.mutationMode = mutationMode

    def __repr__(self):
        return (
            f"GAConfig(generations={self.generations}, "
            f"tournamentSize={self.tournamentSize}, "
            f"elitismSize={self.elitismSize}, "
            f"mutationProbability={self.mutationProbability}, "
            f"rouletteSelection={self.rouletteSelection}, "
            f"crossoverMode={self.crossoverMode}, "
            f"mutationMode={self.mutationMode})"
        )
