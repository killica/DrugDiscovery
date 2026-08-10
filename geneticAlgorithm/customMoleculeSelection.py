"""Lookup molecules from the entire preprocessed molecules catalog."""

from __future__ import annotations

import csv
from pathlib import Path

_DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
FEATURES_META_FILENAME = "features_meta.csv"

_index: dict[str, tuple[str, str]] | None = None


def features_meta_path() -> Path:
    return _DATA_DIR / FEATURES_META_FILENAME


def _load_index() -> dict[str, tuple[str, str]]:
    global _index
    if _index is not None:
        return _index

    path = features_meta_path()
    if not path.exists():
        raise FileNotFoundError(f"Dataset metadata not found: {path}")

    loaded: dict[str, tuple[str, str]] = {}
    with open(path, newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            chembl_id = row["Molecule ChEMBL ID"].strip()
            smiles = row["Smiles"].strip()
            if chembl_id and smiles:
                loaded[chembl_id.upper()] = (chembl_id, smiles)

    _index = loaded
    return _index


def lookup_by_chembl_id(chembl_id: str) -> tuple[str, str] | None:
    """Return ``(chembl_id, smiles)`` if *chembl_id* is in the molecules catalog."""
    key = chembl_id.strip().upper()
    if not key:
        return None
    return _load_index().get(key)
