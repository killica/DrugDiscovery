"""Load saved pIC50 models and predict from SMILES."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import joblib
import numpy as np

from potency.featurizer import MorganFeaturizer, default_fingerprint_config_path


SUPPORTED_MODEL_FAMILIES = ("ridge", "random_forest", "lightgbm")


def default_models_dir() -> Path:
    return default_fingerprint_config_path().parents[1] / "models"


def model_path_for_family(family: str) -> Path:
    if family not in SUPPORTED_MODEL_FAMILIES:
        supported = ", ".join(SUPPORTED_MODEL_FAMILIES)
        raise ValueError(f"Unknown model family {family!r}. Supported: {supported}")
    return default_models_dir() / f"{family}_best.joblib"


def load_model_bundle(path: Path | None = None, *, family: str | None = None) -> dict[str, Any]:
    if path is None:
        if family is None:
            raise ValueError("Provide either path or family")
        path = model_path_for_family(family)
    bundle = joblib.load(path)
    required = {"model_family", "model", "best_params", "fingerprint_config", "metrics"}
    missing = required - bundle.keys()
    if missing:
        raise ValueError(f"Invalid model bundle at {path}: missing keys {sorted(missing)}")
    return bundle


class PotencyPredictor:
    """SMILES → predicted pIC50 using a saved sklearn model bundle."""

    def __init__(self, bundle: dict[str, Any], featurizer: MorganFeaturizer | None = None):
        self.bundle = bundle
        self.model_family: str = bundle["model_family"]
        self.model = bundle["model"]
        self.best_params = bundle["best_params"]
        self.metrics = bundle["metrics"]
        fp_config = bundle["fingerprint_config"]
        self.featurizer = featurizer or MorganFeaturizer(fp_config)
        if self.featurizer.config != fp_config:
            raise ValueError("Featurizer fingerprint config does not match model bundle")

    @classmethod
    def from_bundle_path(cls, path: Path) -> "PotencyPredictor":
        return cls(load_model_bundle(path))

    @classmethod
    def from_family(cls, family: str) -> "PotencyPredictor":
        return cls.from_bundle_path(model_path_for_family(family))

    @classmethod
    def from_default(cls) -> "PotencyPredictor":
        """Best test RMSE model in this project: random forest."""
        return cls.from_family("random_forest")

    def predict_smiles(self, smiles: str) -> float | None:
        x = self.featurizer.featurize_smiles_sparse(smiles)
        if x is None:
            return None
        return float(self.model.predict(x)[0])

    def predict_smiles_batch(self, smiles_list: list[str]) -> list[float | None]:
        x, failed_idx = self.featurizer.featurize_smiles_batch(smiles_list)
        failed = set(failed_idx)
        predictions: list[float | None] = [None] * len(smiles_list)

        if x.shape[0] == 0:
            return predictions

        valid_indices = [i for i in range(len(smiles_list)) if i not in failed]
        y_hat = self.model.predict(x)
        for idx, value in zip(valid_indices, y_hat):
            predictions[idx] = float(value)
        return predictions

    def predict_fingerprint(self, x) -> float:
        """Predict from a precomputed fingerprint row (dense vector or sparse matrix)."""
        if hasattr(x, "reshape") and getattr(x, "ndim", 0) == 1:
            x = x.reshape(1, -1)
        return float(self.model.predict(x)[0])
