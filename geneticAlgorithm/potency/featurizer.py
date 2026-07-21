"""Morgan fingerprint featurization for GA potency prediction.

Must match ``data/02_fingerprints.ipynb`` and ``processed/fingerprint_config.json``
so training and inference use identical bit vectors.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy import sparse


def _project_root() -> Path:
    # geneticAlgorithm/potency/featurizer.py → repo root
    return Path(__file__).resolve().parents[2]


def default_fingerprint_config_path() -> Path:
    return _project_root() / "data" / "processed" / "fingerprint_config.json"


def load_fingerprint_config(path: Path | None = None) -> dict[str, Any]:
    config_path = path or default_fingerprint_config_path()
    with open(config_path, encoding="utf-8") as f:
        return json.load(f)


def build_morgan_generator(config: dict[str, Any]):
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=config["radius"],
        fpSize=config["n_bits"],
        includeChirality=config["use_chirality"],
        useBondTypes=config["use_bond_types"],
    )


class MorganFeaturizer:
    """SMILES → Morgan bit vector using a fixed fingerprint config."""

    def __init__(self, config: dict[str, Any]):
        self.config = dict(config)
        self.n_bits = int(config["n_bits"])
        self._generator = build_morgan_generator(config)

    @classmethod
    def from_config_path(cls, path: Path | None = None) -> "MorganFeaturizer":
        return cls(load_fingerprint_config(path))

    @classmethod
    def from_default(cls) -> "MorganFeaturizer":
        return cls.from_config_path()

    def _fp_to_array(self, bit_vect) -> np.ndarray:
        arr = np.zeros((self.n_bits,), dtype=np.float64)
        DataStructs.ConvertToNumpyArray(bit_vect, arr)
        return arr

    def featurize_mol(self, mol) -> np.ndarray | None:
        if mol is None:
            return None
        return self._fp_to_array(self._generator.GetFingerprint(mol))

    def featurize_smiles(self, smiles: str) -> np.ndarray | None:
        mol = Chem.MolFromSmiles(smiles)
        return self.featurize_mol(mol)

    def featurize_smiles_sparse(self, smiles: str) -> sparse.csr_matrix | None:
        row = self.featurize_smiles(smiles)
        if row is None:
            return None
        return sparse.csr_matrix(row.reshape(1, -1))

    def featurize_smiles_batch(self, smiles_list: list[str]) -> tuple[sparse.csr_matrix, list[int]]:
        """Return (X, failed_indices). Rows align with input order; failed rows are omitted."""
        rows = []
        failed_idx = []
        for i, smiles in enumerate(smiles_list):
            row = self.featurize_smiles(smiles)
            if row is None:
                failed_idx.append(i)
                continue
            rows.append(row)
        if not rows:
            return sparse.csr_matrix((0, self.n_bits)), failed_idx
        return sparse.csr_matrix(np.vstack(rows)), failed_idx
