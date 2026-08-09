"""Paths and I/O for the GA molecule catalog.

- ``molecules_catalog.json`` — fixed seed set (generated from the training dataset)
- ``molecules_user.json`` — molecules added by the user during a session
"""

from __future__ import annotations

import json
from pathlib import Path

_DATA_DIR = Path(__file__).resolve().parent.parent / "data"
CATALOG_SEED_FILENAME = "molecules_catalog.json"
CATALOG_USER_FILENAME = "molecules_user.json"


def molecules_catalog_path() -> Path:
    return _DATA_DIR / CATALOG_SEED_FILENAME


def molecules_user_path() -> Path:
    return _DATA_DIR / CATALOG_USER_FILENAME


def _load_json_list(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as file:
        return json.load(file)


def _write_json_list(path: Path, data: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as file:
        json.dump(data, file, indent=4)
        file.write("\n")


def load_catalog() -> list[dict]:
    """Seed catalog first, then user-added molecules."""
    return _load_json_list(molecules_catalog_path()) + _load_json_list(molecules_user_path())


def append_to_catalog(smiles: str, description: str) -> None:
    """Append a user-added molecule (does not modify the seed catalog)."""
    path = molecules_user_path()
    data = _load_json_list(path)
    data.append({"SMILES": smiles, "Description": description})
    _write_json_list(path, data)
