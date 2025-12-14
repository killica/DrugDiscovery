class GAConfig:
    def __init__(
        self,
        generations: int = 100,
        tournamentSize: int = 4,
        elitismSize: int = 1,
        mutationProbability: float = 0.05,
        rouletteSelection: bool = False
    ):
        self.generations = generations
        self.tournamentSize = tournamentSize
        self.elitismSize = elitismSize
        self.mutationProbability = mutationProbability
        self.rouletteSelection = rouletteSelection

    def __repr__(self):
        return (
            f"GAConfig(generations={self.generations}, "
            f"tournamentSize={self.tournamentSize}, "
            f"elitismSize={self.elitismSize}, "
            f"mutationProbability={self.mutationProbability}, "
            f"rouletteSelection={self.rouletteSelection})"
        )
