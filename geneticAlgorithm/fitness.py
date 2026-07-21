from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import QED

from potency.predictor import PotencyPredictor

# Must match fitness radio ids in main.py (0 = QED, 1–3 = saved models).
FITNESS_MODE_QED = 0
FITNESS_MODE_RANDOM_FOREST = 1
FITNESS_MODE_LIGHTGBM = 2
FITNESS_MODE_RIDGE = 3

MODE_TO_MODEL_FAMILY = {
    FITNESS_MODE_RANDOM_FOREST: "random_forest",
    FITNESS_MODE_LIGHTGBM: "lightgbm",
    FITNESS_MODE_RIDGE: "ridge",
}

MODE_LABELS = {
    FITNESS_MODE_QED: "QED",
    FITNESS_MODE_RANDOM_FOREST: "Random Forest (pIC50)",
    FITNESS_MODE_LIGHTGBM: "LightGBM (pIC50)",
    FITNESS_MODE_RIDGE: "Ridge (pIC50)",
}

INVALID_FITNESS = 0.0


def format_fitness_display(mode: int, value: float) -> str:
    if mode == FITNESS_MODE_QED:
        return f"QED: {value:.4f}"
    return f"pIC50: {value:.4f}"

_predictor_cache: dict[str, PotencyPredictor] = {}


def is_potency_mode(mode: int) -> bool:
    return mode in MODE_TO_MODEL_FAMILY


def mode_label(mode: int) -> str:
    try:
        return MODE_LABELS[mode]
    except KeyError as exc:
        raise ValueError(f"Unknown fitness mode: {mode}") from exc


def get_potency_predictor(mode: int) -> PotencyPredictor:
    family = MODE_TO_MODEL_FAMILY[mode]
    if family not in _predictor_cache:
        _predictor_cache[family] = PotencyPredictor.from_family(family)
    return _predictor_cache[family]


def clear_predictor_cache() -> None:
    """Drop cached models (useful for tests or reloading after retraining)."""
    _predictor_cache.clear()


def _mol_to_smiles(molecule) -> str | None:
    if molecule is None:
        return None
    try:
        return Chem.MolToSmiles(molecule)
    except Exception:
        return None


class Fitness:
    def __init__(
        self,
        molecule,
        weights=(0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95),
        mode: int = FITNESS_MODE_QED,
        smiles: str | None = None,
        predictor: PotencyPredictor | None = None,
    ):
        self.molecule = molecule
        self.weights = weights
        self.mode = mode
        self.smiles = smiles if smiles is not None else _mol_to_smiles(molecule)
        self._predictor = predictor

    def qed(self) -> float:
        if self.molecule is None:
            return INVALID_FITNESS
        return QED.qed(self.molecule, self.weights)

    def potency(self) -> float:
        if not is_potency_mode(self.mode):
            raise ValueError(
                f"potency() requires a model fitness mode, got {mode_label(self.mode)}"
            )
        if not self.smiles:
            return INVALID_FITNESS
        predictor = self._predictor or get_potency_predictor(self.mode)
        prediction = predictor.predict_smiles(self.smiles)
        return prediction if prediction is not None else INVALID_FITNESS

    def value(self) -> float:
        """Active fitness score for the configured mode (higher is better)."""
        if self.mode == FITNESS_MODE_QED:
            return self.qed()
        return self.potency()

    def score(self) -> float:
        return self.value()
