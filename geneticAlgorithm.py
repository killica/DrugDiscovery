import random
import re
import sys
import os
import contextlib
from rdkit import Chem
import mutationInfo

def geneticAlgorithm(population, onlyOneGeneration, numberOfGenerations, rouletteSelection, tournamentSize, elitismSize, mutationProbability, mi, individualLabel, individualProgress):
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

            mutation(newPopulation[j], mutationProbability, mi)
            mutation(newPopulation[j+1], mutationProbability, mi)

            newPopulation[j].setDescription("")
            newPopulation[j].calcFitness()

            newPopulation[j+1].setDescription("")
            newPopulation[j+1].calcFitness()

            individualLabel.setText(f"Individual: {j+1}/{len(population)}")
            individualProgress.setValue(j+1)

        population = newPopulation.copy()
        individualLabel.setText(f"Individual: {j+2}/{len(population)}")
        individualProgress.setValue(j+2)
    
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

# Function to suppress the RDKit warnings and errors
@contextlib.contextmanager
def suppress_rdkit_warnings():
    # Save the current stdout and stderr
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    # Redirect stdout and stderr to devnull (i.e., ignore them)
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')
    try:
        yield  # Continue execution
    finally:
        # Restore the original stdout and stderr
        sys.stdout = original_stdout
        sys.stderr = original_stderr

def isValidSmiles(smiles):
    # Use context manager to suppress output during RDKit validation
    with suppress_rdkit_warnings():
        mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def crossover(parent1, parent2, child1, child2):
    # Assuming smiles strings have more than 1 character (in order to be eligible for crossover)
    smiles1 = parent1.getSmiles()
    smiles2 = parent2.getSmiles()
    MAXITERS = 2000
    i = 0
    while i < MAXITERS:
        i += 1
        cp1 = random.randrange(1, len(smiles1))
        cp2 = random.randrange(1, len(smiles2))

        childSmiles1 = smiles1[:cp1] + smiles2[cp2:]
        childSmiles2 = smiles1[cp1:] + smiles2[:cp2]

        if isValidSmiles(childSmiles1) and isValidSmiles(childSmiles2):
            child1.setSmiles(childSmiles1)
            child2.setSmiles(childSmiles2)
            break
    
    if i == MAXITERS:
        print('Maximum number of iterations for crossover of following molecules exceeded:\n')
        print(f'{smiles1}\n{smiles2}\n')
        print('Passing them to the new generation.\n')
        child1.setSmiles(smiles1)
        child2.setSmiles(smiles2)

def mutation(individual, mutationProbability, mi):
    if random.random() > mutationProbability:
        return

    # Mutation will take place
    mutationType = random.randrange(0, 2)
    # 0 - atom switch
    # 1 - group switch
    # 2 - insertion of an atom or a group - future idea
    # 3 - deletion of an atom or a group - future idea

    if mutationType == 0:
        atomSwitchMutation(individual, mi)
    elif mutationType == 1:
        groupSwitchMutation(individual, mi)
    # elif mutationType == 2:
    #     insertionMutation(individual, mi)
    # else:
    #     deletionMutation(individual, mi)

    
def atomSwitchMutation(individual, mi):
    smiles = individual.getSmiles()
    heteroAtoms = ['O', 'S', 'N', 'P']
    indicesOfHeteroAtoms = []
    for i, ch in enumerate(smiles):
        if ch in heteroAtoms:
            indicesOfHeteroAtoms.append(i)

    if len(indicesOfHeteroAtoms) == 0:
        # Can't perform hetero atom switching
        groupSwitchMutation(individual, mi)
        return
    randomHeteroIndex = random.randrange(len(indicesOfHeteroAtoms))
    heteroAtom = smiles[indicesOfHeteroAtoms[randomHeteroIndex]]

    rnd = random.random()
    cum_prob = 0 # cumulative probability, similar to roulette selection
    changeWith = heteroAtom
    for (second, prob) in mi.atomSwitchMap[heteroAtom]:
        cum_prob += prob
        if rnd <= cum_prob:
            changeWith = second
            break
    smiles = smiles[:indicesOfHeteroAtoms[randomHeteroIndex]] + changeWith + smiles[indicesOfHeteroAtoms[randomHeteroIndex] + 1:]
    individual.setSmiles(smiles)

    with open('log.txt', 'a') as file:
        file.write(f"Changing:{heteroAtom} with {changeWith}, new smiles: {smiles}")

def groupSwitchMutation(individual, mi):
    smiles = individual.getSmiles()
    modifiedSmiles = []

    # try with all known group mutations, and pick a random one
    for replaceFrom, options in mi.groupSwitchMap.items():
        for replaceWith in options:
            # Find all the start positions of replaceFrom in smiles
            matches = [match.start() for match in re.finditer(re.escape(replaceFrom), smiles)]
            if not matches:
                continue
            
            # Generate all replacements by replacing each occurrence of smilesFrom
            for match in matches:
                # Replace only the specific occurrence at the current position
                modSmiles = smiles[:match] + replaceWith + smiles[match + len(replaceFrom):]
                modifiedSmiles.append(modSmiles)

    newSmiles = smiles
    if len(modifiedSmiles) > 0:
        newSmiles = random.choice(modifiedSmiles)
    individual.setSmiles(newSmiles)

    with open('log.txt', 'a') as file:
        file.write(f"Changing:{smiles} with {newSmiles}\n")
    
def insertionMutation(individual, mi):
    pass

def deletionMutation(individual, mi):
    pass


