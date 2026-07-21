import random
import re
import sys
import os
import copy
import contextlib
from rdkit import Chem
from rdkit.Chem import BRICS
import selfies as sf
import mutationInfo
from GAConfig import CrossoverMode, MutationMode
from PyQt5.QtWidgets import QApplication
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

_crossover_stats = {
    mode: {"attempts": 0, "successes": 0, "failures": 0}
    for mode in CrossoverMode
}
_crossover_gen_stats = {
    mode: {"attempts": 0, "successes": 0, "failures": 0}
    for mode in CrossoverMode
}
_crossover_run_generation = 0


def _empty_mutation_stats():
    return {
        "invocations": 0,
        "successes": 0,
        "fallbacks": 0,
        "skipped_invalid": 0,
        "unchanged": 0,
        "total_inner_attempts": 0,
        "inner_attempts_on_success": 0,
        "by_kind": {},
    }


_mutation_stats = {
    mode: _empty_mutation_stats() for mode in MutationMode
}
_mutation_gen_stats = {
    mode: _empty_mutation_stats() for mode in MutationMode
}


def reset_mutation_gen_stats():
    for mode in MutationMode:
        _mutation_gen_stats[mode] = _empty_mutation_stats()


def reset_mutation_stats():
    reset_mutation_gen_stats()
    reset_brics_fragment_pool()
    for mode in MutationMode:
        _mutation_stats[mode] = _empty_mutation_stats()


def _mutation_label(mode: MutationMode) -> str:
    return f"mutation_{mode.name.lower()}"


def _mutation_stats_targets(mode: MutationMode):
    return _mutation_stats[mode], _mutation_gen_stats[mode]


def _mutation_record_invocation(mode: MutationMode):
    for stats in _mutation_stats_targets(mode):
        stats["invocations"] += 1


def _mutation_record_success(mode: MutationMode, inner_attempts: int, kind=None):
    for stats in _mutation_stats_targets(mode):
        stats["successes"] += 1
        stats["total_inner_attempts"] += inner_attempts
        stats["inner_attempts_on_success"] += inner_attempts
        if kind is not None:
            stats["by_kind"][kind] = stats["by_kind"].get(kind, 0) + 1


def _mutation_record_fallback(mode: MutationMode, inner_attempts: int):
    for stats in _mutation_stats_targets(mode):
        stats["fallbacks"] += 1
        stats["total_inner_attempts"] += inner_attempts


def _mutation_record_skipped(mode: MutationMode):
    for stats in _mutation_stats_targets(mode):
        stats["skipped_invalid"] += 1


def _mutation_record_unchanged(mode: MutationMode):
    for stats in _mutation_stats_targets(mode):
        stats["unchanged"] += 1


def _print_mutation_stats(mode: MutationMode):
    stats = _mutation_stats[mode]
    invocations = stats["invocations"]
    if invocations == 0:
        return
    successes = stats["successes"]
    rate = 100.0 * successes / invocations
    avg_attempts = stats["inner_attempts_on_success"] / successes if successes else 0.0
    print(
        f"[{_mutation_label(mode)}] run summary: {successes}/{invocations} successes "
        f"({rate:.1f}%), {stats['fallbacks']} fallbacks, "
        f"{stats['skipped_invalid']} skipped invalid, {stats['unchanged']} unchanged; "
        f"avg attempts on success: {avg_attempts:.1f}"
    )


def _print_mutation_generation_summary(mode: MutationMode):
    stats = _mutation_gen_stats[mode]
    invocations = stats["invocations"]
    label = _mutation_label(mode)
    gen_num = _crossover_run_generation
    if invocations == 0:
        print(f"[{label}] generation {gen_num}: no mutation invocations")
        return
    successes = stats["successes"]
    success_rate = 100.0 * successes / invocations
    fallback_rate = 100.0 * stats["fallbacks"] / invocations
    avg_attempts = stats["inner_attempts_on_success"] / successes if successes else 0.0
    kind_parts = [f"{kind}={count}" for kind, count in sorted(stats["by_kind"].items())]
    kind_str = f", by kind: {', '.join(kind_parts)}" if kind_parts else ""
    print(
        f"[{label}] generation {gen_num} summary: "
        f"{successes}/{invocations} successes ({success_rate:.1f}%), "
        f"{stats['fallbacks']} fallbacks ({fallback_rate:.1f}%), "
        f"{stats['skipped_invalid']} skipped invalid, {stats['unchanged']} unchanged; "
        f"avg attempts on success: {avg_attempts:.1f}{kind_str}"
    )


def get_crossover_generation_stats(mode: CrossoverMode):
    return copy.deepcopy(_crossover_gen_stats[mode])


def get_mutation_generation_stats(mode: MutationMode):
    return copy.deepcopy(_mutation_gen_stats[mode])


def get_crossover_run_stats(mode: CrossoverMode):
    return copy.deepcopy(_crossover_stats[mode])


def get_mutation_run_stats(mode: MutationMode):
    return copy.deepcopy(_mutation_stats[mode])


def reset_crossover_gen_stats():
    for mode in CrossoverMode:
        _crossover_gen_stats[mode] = {"attempts": 0, "successes": 0, "failures": 0}


def reset_crossover_stats():
    global _crossover_run_generation
    _crossover_run_generation = 0
    reset_crossover_gen_stats()
    for mode in CrossoverMode:
        _crossover_stats[mode] = {"attempts": 0, "successes": 0, "failures": 0}


def _crossover_label(mode: CrossoverMode) -> str:
    return f"crossover_{mode.name.lower()}"


def _print_crossover_stats(mode: CrossoverMode):
    stats = _crossover_stats[mode]
    attempts = stats["attempts"]
    if attempts == 0:
        return
    successes = stats["successes"]
    rate = 100.0 * successes / attempts
    print(
        f"[{_crossover_label(mode)}] run success rate: {successes}/{attempts} "
        f"({rate:.1f}%)"
    )


def _print_crossover_generation_summary(mode: CrossoverMode):
    global _crossover_run_generation
    _crossover_run_generation += 1
    gen_num = _crossover_run_generation
    stats = _crossover_gen_stats[mode]
    attempts = stats["attempts"]
    label = _crossover_label(mode)
    if attempts == 0:
        print(f"[{label}] generation {gen_num}: no crossover attempts")
        return
    successes = stats["successes"]
    rate = 100.0 * successes / attempts
    print(f"[{label}] generation {gen_num} summary: {successes}/{attempts} ({rate:.1f}%)")


def _crossover_record_success(mode: CrossoverMode):
    for stats in (_crossover_stats[mode], _crossover_gen_stats[mode]):
        stats["attempts"] += 1
        stats["successes"] += 1


def _crossover_use_parents(mode, reason, smiles1, smiles2, child1, child2):
    for stats in (_crossover_stats[mode], _crossover_gen_stats[mode]):
        stats["attempts"] += 1
        stats["failures"] += 1
    print(f"[{_crossover_label(mode)}] FAILED ({reason}) — using parent structures:")
    print(f"  parent1: {smiles1}")
    print(f"  parent2: {smiles2}")
    _print_crossover_stats(mode)
    child1.setSmiles(smiles1)
    child2.setSmiles(smiles2)

def _process_gui_events():
    """Allow progress labels/bars to repaint while GA runs on the main thread."""
    app = QApplication.instance()
    if app is not None:
        app.processEvents()


def _set_individual_progress(individualLabel, individualProgress, current, total):
    """Update progress widgets; ignore if Qt has already destroyed them."""
    if individualLabel is None or individualProgress is None:
        return
    try:
        total = max(total, 1)
        current = max(0, min(current, total))
        individualLabel.setText(f"Individual: {current}/{total}")
        individualProgress.setMaximum(total)
        individualProgress.setValue(current)
    except RuntimeError:
        pass


def _crossover_for_mode(mode: CrossoverMode):
    if mode == CrossoverMode.GRAPH:
        return crossover_graph
    if mode == CrossoverMode.SELFIES:
        return crossover_selfies
    if mode == CrossoverMode.BRICS:
        return crossover_brics
    return crossover_smiles


def geneticAlgorithm(
    population,
    onlyOneGeneration,
    generations,
    rouletteSelection,
    tournamentSize,
    elitismSize,
    mutationProbability,
    crossoverMode,
    mutationMode,
    mi,
    individualLabel,
    individualProgress,
    cancel_check=None,
    on_generation_start=None,
    on_new_individual=None,
    evolution_stats=None,
    record_initial=False,
):
    crossover = _crossover_for_mode(crossoverMode)
    populationSize = len(population)
    if cancel_check and cancel_check():
        return population
    _set_individual_progress(individualLabel, individualProgress, 0, populationSize)
    _process_gui_events()
    newPopulation = population.copy()
    if onlyOneGeneration:
        generations = 1
    
    if elitismSize % 2 != populationSize % 2:
        elitismSize += 1

    if mutationMode == MutationMode.BRICS:
        init_brics_fragment_pool(population)

    if evolution_stats is not None and record_initial:
        evolution_stats.record_generation(
            population,
            crossover_mode=crossoverMode,
            mutation_mode=mutationMode,
            crossover_stats=get_crossover_generation_stats(crossoverMode),
            mutation_stats=get_mutation_generation_stats(mutationMode),
        )

    for gen_idx in range(generations):
        if cancel_check and cancel_check():
            return population
        reset_crossover_gen_stats()
        reset_mutation_gen_stats()
        # current population is already sorted
        newPopulation[:elitismSize] = population[:elitismSize]
        if on_generation_start is not None:
            on_generation_start(gen_idx, populationSize)
        tmp = 0
        _process_gui_events()

        for i in range(elitismSize):
            if cancel_check and cancel_check():
                return population
            if on_new_individual is not None:
                on_new_individual(newPopulation[i], i, gen_idx)
            _process_gui_events()

        for j in range(elitismSize, populationSize, 2):
            if cancel_check and cancel_check():
                return population
            parent1 = selection(population, rouletteSelection, tournamentSize)
            parent2 = selection(population, rouletteSelection, tournamentSize) # TODO: Parents be different

            crossover(parent1, parent2, child1=newPopulation[j], child2=newPopulation[j + 1])

            mutation(newPopulation[j], mutationProbability, mutationMode, mi)
            mutation(newPopulation[j + 1], mutationProbability, mutationMode, mi)

            newPopulation[j].setDescription("")
            newPopulation[j].calcFitness()

            newPopulation[j+1].setDescription("")
            newPopulation[j+1].calcFitness()

            _set_individual_progress(individualLabel, individualProgress, j + 1, len(population))
            if on_new_individual is not None:
                on_new_individual(newPopulation[j], j, gen_idx)
            _process_gui_events()
            if on_new_individual is not None:
                on_new_individual(newPopulation[j + 1], j + 1, gen_idx)
            _process_gui_events()

            tmp = j

        population = newPopulation.copy()
        if evolution_stats is not None:
            evolution_stats.record_generation(
                population,
                crossover_mode=crossoverMode,
                mutation_mode=mutationMode,
                crossover_stats=get_crossover_generation_stats(crossoverMode),
                mutation_stats=get_mutation_generation_stats(mutationMode),
            )
        _set_individual_progress(individualLabel, individualProgress, tmp + 2, len(population))
        _process_gui_events()
        _print_crossover_generation_summary(crossoverMode)
        _print_mutation_generation_summary(mutationMode)

    if evolution_stats is not None:
        evolution_stats.update_operators_summary(
            crossoverMode,
            mutationMode,
            get_crossover_run_stats(crossoverMode),
            get_mutation_run_stats(mutationMode),
        )
    _print_mutation_stats(mutationMode)
    return population

def selection(population, rouletteSelection, tournamentSize):
    if rouletteSelection:
        return rouletteWheelSelection(population)
    return tournamentSelection(population, tournamentSize)

def rouletteWheelSelection(population):
    totalFitness = sum(individual.getFitness() for individual in population)
    probabilitiesPartialSums = []
    tmp = 0
    for individual in population:
        tmp += individual.getFitness() / totalFitness
        probabilitiesPartialSums.append(tmp)
    
    randomValue = random.random()
    
    for i, partialSum in enumerate(probabilitiesPartialSums):
        if randomValue <= partialSum:
            return population[i]
    
    return population[-1]

def tournamentSelection(population, tournamentSize):
    populationSample = random.sample(population, k = tournamentSize)
    bestFitness = -1
    bestIndividual = None
    for individual in populationSample:
        if individual.getFitness() > bestFitness:
            bestFitness = individual.getFitness()
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


def _brics_build_random(frags, max_candidates=32):
    """Return one BRICS-rebuilt molecule without enumerating all combinations."""
    candidates = []
    for mol in BRICS.BRICSBuild(frags):
        candidates.append(mol)
        if len(candidates) >= max_candidates:
            break
    if not candidates:
        return None
    return random.choice(candidates)


_DEFAULT_BRICS_SEED_SMILES = [
    "c1ccccc1",
    "CCO",
    "CC(=O)O",
    "c1ccncc1",
    "CC(C)C",
    "CCN(CC)CC",
    "c1ccc(O)cc1",
]

_brics_fragment_pool = []


def build_brics_fragment_pool(smiles_list):
    pool = set()
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        try:
            pool.update(BRICS.BRICSDecompose(mol))
        except Exception:
            pass
    return list(pool)


def _default_brics_fragment_pool():
    return build_brics_fragment_pool(_DEFAULT_BRICS_SEED_SMILES)


def init_brics_fragment_pool(population):
    global _brics_fragment_pool
    smiles_list = [individual.getSmiles() for individual in population]
    pool = set(build_brics_fragment_pool(smiles_list))
    pool.update(_default_brics_fragment_pool())
    _brics_fragment_pool = list(pool)
    print(f"[mutation_brics] fragment pool size: {len(_brics_fragment_pool)}")


def reset_brics_fragment_pool():
    global _brics_fragment_pool
    _brics_fragment_pool = []


def crossover_smiles(parent1, parent2, child1, child2, MAX_ITERS = 2000):
    mode = CrossoverMode.SMILES
    smiles1 = parent1.getSmiles()
    smiles2 = parent2.getSmiles()
    n1 = len(smiles1)
    n2 = len(smiles2)

    if n1 < 2 or n2 < 2:
        _crossover_use_parents(
            mode, "molecule too short for single-point crossover", smiles1, smiles2, child1, child2
        )
        return

    i = 0
    while i < MAX_ITERS:
        i += 1
        cp1 = random.randrange(1, n1)
        cp2 = random.randrange(1, n2)

        childSmiles1 = smiles1[:cp1] + smiles2[cp2:]
        childSmiles2 = smiles1[cp1:] + smiles2[:cp2]

        if isValidSmiles(childSmiles1) and isValidSmiles(childSmiles2):
            child1.setSmiles(childSmiles1)
            child2.setSmiles(childSmiles2)
            _crossover_record_success(mode)
            return

    if n1 < 3 or n2 < 3:
        _crossover_use_parents(
            mode,
            f"single-point failed after {MAX_ITERS} attempts; molecule too short for two-point crossover",
            smiles1,
            smiles2,
            child1,
            child2,
        )
        return

    i = 0
    while i < MAX_ITERS:
        i += 1
        first1 = random.randrange(1, n1)
        if first1 < n1 - 1:
            first2 = random.randrange(first1 + 1, n1)
        else:
            first2 = random.randrange(1, first1)
            # first2 < first1, so swap them
            first1, first2 = first2, first1

        second1 = random.randrange(1, n2)
        if second1 < n2 - 1:
            second2 = random.randrange(second1 + 1, n2)
        else:
            second2 = random.randrange(1, second1)
            # second2 < second1, so swap them
            second1, second2 = second2, second1

        childSmiles1 = smiles1[:first1] + smiles2[second1:second2] + smiles1[first2:]
        childSmiles2 = smiles2[:second1] + smiles1[first1:first2] + smiles2[second2:]

        if isValidSmiles(childSmiles1) and isValidSmiles(childSmiles2):
            child1.setSmiles(childSmiles1)
            child2.setSmiles(childSmiles2)
            _crossover_record_success(mode)
            return

    _crossover_use_parents(
        mode,
        f"no valid recombination after {MAX_ITERS} single-point and {MAX_ITERS} two-point attempts",
        smiles1,
        smiles2,
        child1,
        child2,
    )

def crossover_selfies(parent1, parent2, child1, child2, MAX_ITERS = 20):
    mode = CrossoverMode.SELFIES
    smiles1 = parent1.getSmiles()
    smiles2 = parent2.getSmiles()
    try:
        selfies1 = sf.encoder(smiles1)
        selfies2 = sf.encoder(smiles2)

    except Exception:
        _crossover_use_parents(mode, "SELFIES encoding failed", smiles1, smiles2, child1, child2)
        return

    tokens1 = list(sf.split_selfies(selfies1))
    tokens2 = list(sf.split_selfies(selfies2))

    if len(tokens1) < 2 or len(tokens2) < 2:
        _crossover_use_parents(mode, "molecule too small for crossover", smiles1, smiles2, child1, child2)
        return

    for _ in range(MAX_ITERS):

        # Select random crossover points
        cut1 = random.randint(1, len(tokens1) - 1)
        cut2 = random.randint(1, len(tokens2) - 1)

        # Create children token sequences
        child1_tokens = tokens1[:cut1] + tokens2[cut2:]
        child2_tokens = tokens2[:cut2] + tokens1[cut1:]

        child1_selfies = "".join(child1_tokens)
        child2_selfies = "".join(child2_tokens)

        try:
            # Decode SELFIES -> SMILES
            child1_smiles = sf.decoder(child1_selfies)
            child2_smiles = sf.decoder(child2_selfies)

            # Validate with RDKit
            mol1 = Chem.MolFromSmiles(child1_smiles)
            mol2 = Chem.MolFromSmiles(child2_smiles)

            if mol1 is None or mol2 is None:
                continue

            # Canonicalize SMILES
            child1_smiles = Chem.MolToSmiles(mol1, canonical=True)
            child2_smiles = Chem.MolToSmiles(mol2, canonical=True)

            child1.setSmiles(child1_smiles)
            child2.setSmiles(child2_smiles)
            _crossover_record_success(mode)
            return

        except Exception:
            continue

    _crossover_use_parents(
        mode,
        f"no valid recombination after {MAX_ITERS} attempts",
        smiles1,
        smiles2,
        child1,
        child2,
    )

def crossover_graph(parent1, parent2, child1, child2, MAX_ITERS = 100):
    mode = CrossoverMode.GRAPH
    smiles1 = parent1.getSmiles()
    smiles2 = parent2.getSmiles()
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        _crossover_use_parents(mode, "invalid parent SMILES", smiles1, smiles2, child1, child2)
        return

    # Get non-ring bonds
    bonds1 = [
        bond.GetIdx()
        for bond in mol1.GetBonds()
        if not bond.IsInRing()
    ]

    bonds2 = [
        bond.GetIdx()
        for bond in mol2.GetBonds()
        if not bond.IsInRing()
    ]

    if not bonds1 or not bonds2:
        _crossover_use_parents(mode, "no cuttable non-ring bonds", smiles1, smiles2, child1, child2)
        return

    for _ in range(MAX_ITERS):
        try:
            # Select random bonds for cutting
            bond_idx1 = random.choice(bonds1)
            bond_idx2 = random.choice(bonds2)

            # Fragment molecules -> add dummy atoms at the cut bonds
            frag_mol1 = Chem.FragmentOnBonds(
                mol1,
                [bond_idx1],
                addDummies=True
            )

            frag_mol2 = Chem.FragmentOnBonds(
                mol2,
                [bond_idx2],
                addDummies=True
            )

            # Extract fragments
            frags1 = Chem.GetMolFrags(
                frag_mol1,
                asMols=True,
                sanitizeFrags=True
            )

            frags2 = Chem.GetMolFrags(
                frag_mol2,
                asMols=True,
                sanitizeFrags=True
            )

            # If cut bonds were not parts of the ring, we expect 2 fragments
            if len(frags1) < 2 or len(frags2) < 2:
                continue

            # Randomly recombine fragments - we group two disconnected fragments into one molecule
            child1_base = Chem.CombineMols(
                frags1[0],
                frags2[1]
            )

            child2_base = Chem.CombineMols(
                frags2[0],
                frags1[1]
            )

            # Editable molecules
            rw_child1 = Chem.RWMol(child1_base)
            rw_child2 = Chem.RWMol(child2_base)

            # Find dummy atoms
            dummy_atoms1 = [
                atom.GetIdx()
                for atom in rw_child1.GetAtoms()
                if atom.GetAtomicNum() == 0
            ]

            dummy_atoms2 = [
                atom.GetIdx()
                for atom in rw_child2.GetAtoms()
                if atom.GetAtomicNum() == 0
            ]

            # Need exactly 2 dummy atoms to reconnect
            if len(dummy_atoms1) != 2 or len(dummy_atoms2) != 2:
                continue

            # Neighbors of dummy atoms - they will be used to reconnect the fragments
            nbrs1 = [
                rw_child1.GetAtomWithIdx(idx).GetNeighbors()[0].GetIdx()
                for idx in dummy_atoms1
            ]

            nbrs2 = [
                rw_child2.GetAtomWithIdx(idx).GetNeighbors()[0].GetIdx()
                for idx in dummy_atoms2
            ]

            # Remove dummy atoms (reverse order!)
            for idx in sorted(dummy_atoms1, reverse=True):
                rw_child1.RemoveAtom(idx)

            for idx in sorted(dummy_atoms2, reverse=True):
                rw_child2.RemoveAtom(idx)

            # Adjust indices after deletion
            def adjust_index(idx, removed):
                shift = sum(1 for r in removed if r < idx)
                return idx - shift

            removed1 = sorted(dummy_atoms1)
            removed2 = sorted(dummy_atoms2)

            nbrs1 = [
                adjust_index(i, removed1)
                for i in nbrs1
            ]

            nbrs2 = [
                adjust_index(i, removed2)
                for i in nbrs2
            ]

            # Reconnect fragments - Single bond is the safest option
            rw_child1.AddBond(
                nbrs1[0],
                nbrs1[1],
                Chem.BondType.SINGLE
            )

            rw_child2.AddBond(
                nbrs2[0],
                nbrs2[1],
                Chem.BondType.SINGLE
            )

            # Final molecules
            child1_mol = rw_child1.GetMol()
            child2_mol = rw_child2.GetMol()

            # Sanitize
            Chem.SanitizeMol(child1_mol)
            Chem.SanitizeMol(child2_mol)

            # Convert to canonical SMILES
            child1_smiles = Chem.MolToSmiles(
                child1_mol,
                canonical=True
            )

            child2_smiles = Chem.MolToSmiles(
                child2_mol,
                canonical=True
            )

            child1.setSmiles(child1_smiles)
            child2.setSmiles(child2_smiles)
            _crossover_record_success(mode)
            return

        except Exception:
            continue

    _crossover_use_parents(
        mode,
        f"no valid recombination after {MAX_ITERS} attempts",
        smiles1,
        smiles2,
        child1,
        child2,
    )

def crossover_brics(parent1, parent2, child1, child2, MAX_ITERS=100, MAX_BRICS_CANDIDATES=10):
    mode = CrossoverMode.BRICS
    smiles1 = parent1.getSmiles()
    smiles2 = parent2.getSmiles()
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        _crossover_use_parents(mode, "invalid parent SMILES", smiles1, smiles2, child1, child2)
        return

    frags1 = list(BRICS.BRICSDecompose(mol1))
    frags2 = list(BRICS.BRICSDecompose(mol2))

    if len(frags1) < 2 or len(frags2) < 2:
        _crossover_use_parents(
            mode, "insufficient BRICS fragments", smiles1, smiles2, child1, child2
        )
        return

    for _ in range(MAX_ITERS):
        _process_gui_events()
        try:
            cut1 = random.randint(1, len(frags1) - 1)
            cut2 = random.randint(1, len(frags2) - 1)

            child1_frag_smiles = frags1[:cut1] + frags2[cut2:]
            child2_frag_smiles = frags2[:cut2] + frags1[cut1:]

            child1_frags = []
            child2_frags = []
            valid = True

            for frag in child1_frag_smiles:
                frag_mol = Chem.MolFromSmiles(frag)
                if frag_mol is None:
                    valid = False
                    break
                child1_frags.append(frag_mol)

            if valid:
                for frag in child2_frag_smiles:
                    frag_mol = Chem.MolFromSmiles(frag)
                    if frag_mol is None:
                        valid = False
                        break
                    child2_frags.append(frag_mol)

            if not valid:
                continue

            child1_mol = _brics_build_random(child1_frags, MAX_BRICS_CANDIDATES)
            child2_mol = _brics_build_random(child2_frags, MAX_BRICS_CANDIDATES)

            if child1_mol is None or child2_mol is None:
                continue

            Chem.SanitizeMol(child1_mol)
            Chem.SanitizeMol(child2_mol)

            child1.setSmiles(Chem.MolToSmiles(child1_mol, canonical=True))
            child2.setSmiles(Chem.MolToSmiles(child2_mol, canonical=True))
            _crossover_record_success(mode)
            return

        except Exception:
            continue

    _crossover_use_parents(
        mode,
        f"no valid recombination after {MAX_ITERS} attempts",
        smiles1,
        smiles2,
        child1,
        child2,
    )

def _apply_smiles_to_individual(individual, smiles):
    if not isValidSmiles(smiles):
        return False
    mol = Chem.MolFromSmiles(smiles)
    individual.setSmiles(Chem.MolToSmiles(mol, canonical=True))
    return True


def _mutate_for_mode(mode: MutationMode):
    if mode == MutationMode.SELFIES:
        return mutate_selfies
    if mode == MutationMode.GRAPH:
        return mutate_graph
    if mode == MutationMode.BRICS:
        return mutate_brics
    return mutate_smiles


def mutation(individual, mutationProbability, mutationMode, mi):
    if random.random() > mutationProbability:
        return
    _mutation_record_invocation(mutationMode)
    mutate = _mutate_for_mode(mutationMode)
    mutate(individual, mi)


def mutate_smiles(individual, mi):
    mode = MutationMode.SMILES
    smiles_before = individual.getSmiles()
    if not isValidSmiles(smiles_before):
        _mutation_record_skipped(mode)
        return

    mutationType = random.randrange(0, 4)
    # 0 - atom switch
    # 1 - group switch
    # 2 - insertion of an atom or a group
    # 3 - deletion of an atom or a group

    if mutationType == 0:
        atomSwitchMutation(individual, mi)
    elif mutationType == 1:
        groupSwitchMutation(individual, mi)
    elif mutationType == 2:
        insertionMutation(individual, mi)
    else:
        deletionMutation(individual, mi)

    if individual.getSmiles() != smiles_before:
        _mutation_record_success(mode, 1)
    else:
        _mutation_record_unchanged(mode)


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
    if not _apply_smiles_to_individual(individual, smiles):
        groupSwitchMutation(individual, mi)

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
    _apply_smiles_to_individual(individual, newSmiles)
    
def insertionMutation(individual, mi):
    smiles = individual.getSmiles()
    n = len(smiles)
    MAX_INSERTION_ATTEMPTS = 200
    randomInsertion = random.choice(mi.insertions)
    i = 0
    while i < MAX_INSERTION_ATTEMPTS:
        i += 1
        randomInsertionPosition = random.randrange(n)
        newSmiles = smiles[:randomInsertionPosition] + randomInsertion + smiles[randomInsertionPosition:]
        if isValidSmiles(newSmiles):
            _apply_smiles_to_individual(individual, newSmiles)
            return

    # Insertion failed, in order to perform any other mutation, call deletionMutation (for example)
    deletionMutation(individual, mi)

def deletionMutation(individual, mi):
    smiles = individual.getSmiles()
    n = len(smiles)
    MAX_BRANCH_DELETION_ATTEMPTS = 20
    MAX_ITERS = 500
    openBracketsPositions = []
    for i, ch in enumerate(smiles):
        if ch == '(':
            openBracketsPositions.append(i)

    # Branches in the molecule exist
    if len(openBracketsPositions) > 0:
        # Attempt deleting random branches of a molecule
        for _ in range(MAX_BRANCH_DELETION_ATTEMPTS):
            startIndex = random.choice(openBracketsPositions)
            openBrackets = 1
            endIndex = None
            # Finding the corresponding closed bracket
            for i in range(startIndex + 1, n):
                if smiles[i] == '(':
                    openBrackets += 1
                if smiles[i] == ')':
                    openBrackets -= 1
                    if openBrackets == 0:
                        endIndex = i
                        break

            newSmiles = smiles[:startIndex] + smiles[endIndex + 1:]
            if isValidSmiles(newSmiles):
                _apply_smiles_to_individual(individual, newSmiles)
                return

    # Else: zero branches in the molecule, or branches can not be deleted
    i = 0
    while i < MAX_ITERS:
        i += 1
        randomIndex = random.randrange(n)
        # We only want to delete atoms
        if smiles[randomIndex].isalpha():
            newSmiles = smiles[:randomIndex] + smiles[randomIndex + 1:]
            if isValidSmiles(newSmiles):
                _apply_smiles_to_individual(individual, newSmiles)
                break
    # Single atom deletion didn't succeed neither, hence, in order not to leave the molecule unchanged, we call
    # atom switch mutation. In case it also fails, that method will call group mutation, which will not fail,
    # since there are certain mutations that will always succeed.
    if i == MAX_ITERS:
        atomSwitchMutation(individual, mi)


def mutate_selfies(individual, mi):
    mode = MutationMode.SELFIES
    smiles = individual.getSmiles()
    if not isValidSmiles(smiles):
        _mutation_record_skipped(mode)
        return
    # Canonicalize input
    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)
    try:
        selfies_str = sf.encoder(smiles)

    except Exception:
        deletionMutation(individual, mi)
        _mutation_record_fallback(mode, 0)
        return

    alphabet = list(sf.get_semantic_robust_alphabet())
    original_tokens = list(sf.split_selfies(selfies_str))
    if len(original_tokens) == 0:
        deletionMutation(individual, mi)
        _mutation_record_fallback(mode, 0)
        return
    
    MAX_MUTATION_ATTEMPTS = 20
    for attempt in range(1, MAX_MUTATION_ATTEMPTS + 1):

        tokens = original_tokens.copy()

        # Choose mutation type (prevent empty molecule by not allowing deletion if there is only one token)
        if len(tokens) > 1:
            mutation_type = random.choice(["replace", "insert", "delete"])
        else:
            mutation_type = random.choice(["replace", "insert"])

        try:
            if mutation_type == "replace":
                idx = random.randrange(len(tokens))
                tokens[idx] = random.choice(alphabet)

            elif mutation_type == "insert":
                idx = random.randrange(len(tokens) + 1)
                tokens.insert(idx, random.choice(alphabet))

            elif mutation_type == "delete":
                idx = random.randrange(len(tokens))
                del tokens[idx]

            # Convert back to SELFIES
            mutated_selfies = "".join(tokens)

            # Decode SELFIES -> SMILES
            mutated_smiles = sf.decoder(mutated_selfies)

            # Validate with RDKit
            mutated_mol = Chem.MolFromSmiles(mutated_smiles)

            if mutated_mol is None:
                continue

            # Sanitize
            Chem.SanitizeMol(mutated_mol)

            # Canonical SMILES
            mutated_smiles = Chem.MolToSmiles(mutated_mol, canonical=True)

            # Avoid identical molecule
            if mutated_smiles == smiles:
                continue

            individual.setSmiles(mutated_smiles)
            _mutation_record_success(mode, attempt, mutation_type)
            return

        except Exception:
            continue

    deletionMutation(individual, mi)
    _mutation_record_fallback(mode, MAX_MUTATION_ATTEMPTS)

_GRAPH_ALLOWED_ATOMS = [
    6,   # C
    7,   # N
    8,   # O
    9,   # F
    16,  # S
    17,  # Cl
]

def mutate_graph(individual, mi):
    mode = MutationMode.GRAPH
    smiles = individual.getSmiles()
    if not isValidSmiles(smiles):
        _mutation_record_skipped(mode)
        return

    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol, canonical=True)

    MAX_MUTATION_ATTEMPTS = 50
    for attempt in range(1, MAX_MUTATION_ATTEMPTS + 1):
        try:
            rwmol = Chem.RWMol(mol)
            mutation_type = random.choice(["atom", "bond"])

            if mutation_type == "atom":
                atoms = [
                    atom.GetIdx()
                    for atom in rwmol.GetAtoms()
                    if atom.GetAtomicNum() > 1
                ]
                if not atoms:
                    continue

                idx = random.choice(atoms)
                rwmol.GetAtomWithIdx(idx).SetAtomicNum(
                    random.choice(_GRAPH_ALLOWED_ATOMS)
                )

            else:
                bonds = list(rwmol.GetBonds())
                if not bonds:
                    continue

                bond_idx = random.choice(bonds).GetIdx()
                new_type = random.choice([
                    Chem.BondType.SINGLE,
                    Chem.BondType.DOUBLE,
                ])
                rwmol.GetBondWithIdx(bond_idx).SetBondType(new_type)

            mutated = rwmol.GetMol()
            # Sanitization fails if the mutated molecule is invalid
            Chem.SanitizeMol(mutated)
            new_smiles = Chem.MolToSmiles(mutated, canonical=True)

            if new_smiles != smiles:
                individual.setSmiles(new_smiles)
                _mutation_record_success(mode, attempt, mutation_type)
                return

        except Exception:
            continue

    deletionMutation(individual, mi)
    _mutation_record_fallback(mode, MAX_MUTATION_ATTEMPTS)


def mutate_brics(individual, mi):
    mode = MutationMode.BRICS
    smiles = individual.getSmiles()
    if not isValidSmiles(smiles):
        _mutation_record_skipped(mode)
        return

    mol = Chem.MolFromSmiles(smiles)
    original_smiles = Chem.MolToSmiles(mol, canonical=True)

    try:
        fragments = list(BRICS.BRICSDecompose(mol))
    except Exception:
        deletionMutation(individual, mi)
        _mutation_record_fallback(mode, 0)
        return

    if len(fragments) == 0:
        deletionMutation(individual, mi)
        _mutation_record_fallback(mode, 0)
        return

    fragment_pool = _brics_fragment_pool
    if len(fragment_pool) == 0:
        deletionMutation(individual, mi)
        _mutation_record_fallback(mode, 0)
        return

    MAX_MUTATION_ATTEMPTS = 20
    MAX_BRICS_CANDIDATES = 32

    for attempt in range(1, MAX_MUTATION_ATTEMPTS + 1):
        try:
            mutated_fragments = fragments.copy()
            idx = random.randrange(len(mutated_fragments))
            replacement = random.choice(fragment_pool)

            if replacement == mutated_fragments[idx]:
                continue

            mutated_fragments[idx] = replacement

            fragment_mols = []
            valid = True
            for frag in mutated_fragments:
                frag_mol = Chem.MolFromSmiles(frag)
                if frag_mol is None:
                    valid = False
                    break
                fragment_mols.append(frag_mol)

            if not valid:
                continue

            candidate = _brics_build_random(fragment_mols, MAX_BRICS_CANDIDATES)
            if candidate is None:
                continue

            Chem.SanitizeMol(candidate)
            candidate_smiles = Chem.MolToSmiles(candidate, canonical=True)

            if candidate_smiles != original_smiles:
                individual.setSmiles(candidate_smiles)
                _mutation_record_success(mode, attempt, "replace")
                return

        except Exception:
            continue

    deletionMutation(individual, mi)
    _mutation_record_fallback(mode, MAX_MUTATION_ATTEMPTS)