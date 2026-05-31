"""Population diversity metrics from pairwise Tanimoto similarity."""

import statistics

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import FingerprintSimilarity

MORGAN_RADIUS = 2
MORGAN_N_BITS = 2048


def pairwise_tanimoto_similarities(smiles_list):
    """Return Tanimoto similarities for all unique molecule pairs."""
    mols = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mols.append(mol)

    if len(mols) < 2:
        return []

    fingerprints = [
        rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol,
            MORGAN_RADIUS,
            nBits=MORGAN_N_BITS,
        )
        for mol in mols
    ]

    similarities = []
    for i in range(len(fingerprints)):
        for j in range(i + 1, len(fingerprints)):
            similarities.append(
                float(FingerprintSimilarity(fingerprints[i], fingerprints[j]))
            )
    return similarities


def summarize_diversity(similarities):
    if not similarities:
        return {
            "pair_count": 0,
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
        }
    return {
        "pair_count": len(similarities),
        "mean": round(sum(similarities) / len(similarities), 4),
        "std": round(statistics.pstdev(similarities), 4),
        "min": round(min(similarities), 4),
        "max": round(max(similarities), 4),
    }


def compute_population_diversity(population):
    smiles_list = [individual.getSmiles() for individual in population]
    similarities = pairwise_tanimoto_similarities(smiles_list)
    return {
        "pairwise_tanimoto": similarities,
        "summary": summarize_diversity(similarities),
    }
