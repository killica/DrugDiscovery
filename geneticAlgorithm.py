import random
from rdkit import Chem

def geneticAlgorithm(population, onlyOneGeneration, numberOfGenerations, rouletteSelection, tournamentSize, elitismSize, mutationProbability):
    populationSize = len(population)
    newPopulation = population.copy()
    if onlyOneGeneration:
        numberOfGenerations = 1
    
    for _ in range(numberOfGenerations):
        # current population is already sorted
        newPopulation[:elitismSize] = population[:elitismSize]
        for j in range(elitismSize, populationSize - 1):
            parent1 = selection(population, rouletteSelection, tournamentSize)
            parent2 = selection(population, rouletteSelection, tournamentSize) # TODO: Parents be different

            crossover(parent1, parent2, child1 = newPopulation[j], child2 = newPopulation[j+1])

            # mutation(newPopulation[j], mutationProbability)
            # mutation(newPopulation[j+1], mutationProbability)

            newPopulation[j].setDescription("")
            newPopulation[j].calcFitness()

            newPopulation[j+1].setDescription("")
            newPopulation[j+1].calcFitness()

        population = newPopulation.copy()
    
    return population

def selection(population, rouletteSelection, tournamentSize):
    if rouletteSelection:
        return rouletteWheelSelection(population)
    return tournamentSelection(population, tournamentSize)

def rouletteWheelSelection(population):
    totalFitness = sum([individual.getQED() for individual in population])
    probabilitiesPartialSums = []
    tmp = 0
    for individual in population:
        tmp += individual.getQED() / totalFitness
        probabilitiesPartialSums.append(tmp)
    
    randomValue = random.random()
    
    for i, partialSum in enumerate(probabilitiesPartialSums):
        if randomValue <= partialSum:
            return population[i]
    
    return population[-1]

def tournamentSelection(population, tournamentSize):
    populationSample = random.sample(population, k = tournamentSize)
    bestQED = -1
    bestIndividual = None
    for individual in populationSample:
        if individual.getQED() > bestQED:
            bestQED = individual.getQED()
            bestIndividual = individual
    return bestIndividual

def crossover(parent1, parent2, child1, child2):
    # Assuming smiles strings have more than 1 character (in order to be eligible for crossover)
    smiles1 = parent1.getSmiles()
    smiles2 = parent2.getSmiles()

    while True:
        cp1 = random.randrange(1, len(smiles1))
        cp2 = random.randrange(1, len(smiles2))

        childSmiles1 = smiles1[:cp1] + smiles2[cp2:]
        childSmiles2 = smiles1[cp1:] + smiles2[:cp2]

        mol1 = Chem.MolFromSmiles(childSmiles1)
        mol2 = Chem.MolFromSmiles(childSmiles2)

        if mol1 is not None and mol2 is not None:
            child1.setSmiles(childSmiles1)
            child2.setSmiles(childSmiles2)
            break



def mutation(individual, mutationProbability):
    pass
    