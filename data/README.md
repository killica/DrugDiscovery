# EGFR pIC50 Activity Prediction

Machine Learning course project, Faculty of Mathematics, University of Belgrade.

Part of the wider [DrugDiscovery](../README.md) repository (Master's thesis + genetic algorithm for in silico drug exploration).

## Author

- Lazar Savić

## Project description

**EGFR** (Epidermal Growth Factor Receptor) is a protein kinase that drives cell growth and a major oncology drug target — inhibitors block its signaling to slow tumor proliferation.

**IC50** is the compound concentration (here in **nM**) needed to inhibit the target by 50%. Lower IC50 → more potent molecule. Databases often use the log-scale **pIC50**:

```
pIC50 = -log₁₀(IC50 in M) = 9 - log₁₀(IC50_nM)
```

We **predict pIC50 from molecular structure** (regression on Morgan fingerprints, explained in [02_fingerprints.ipynb](./02_fingerprints.ipynb)). Higher pIC50 means lower IC50 and thus **greater potency**.

## Directory structure

```
data/
├── README.md
├── Masinsko ucenje prezentacija.pdf   ← project presentation
├── molecules.csv                      ← raw ChEMBL export
├── 01_preprocessing.ipynb             ← data analysis, cleaning, aggregation
├── 02_fingerprints.ipynb              ← Morgan fingerprint generation
├── 03_models.ipynb                    ← training, evaluation, model saving
├── assets/                            ← figures referenced in notebooks
├── processed/                         ← intermediate pipeline artifacts
└── models/                            ← trained models, metrics, tuning logs
```

## Dataset

Bioactivity data for **EGFR** (`CHEMBL203`) from [ChEMBL](https://www.ebi.ac.uk/chembl/), stored in `molecules.csv`. Full exploratory analysis and cleaning are in [01_preprocessing.ipynb](./01_preprocessing.ipynb).

## How to run

Run the notebooks **in order** from the `data/` directory (`01 → 02 → 03`). Paths in the notebooks are relative to `data/` (e.g. `molecules.csv`, `processed/`, `models/`).

## Environment setup

**Prerequisites:** Python **3.9+**, `git`, and a Jupyter environment (JupyterLab, Notebook, or VS Code/Cursor with a Python kernel).

From the **repository root**:

```bash
git clone https://github.com/killica/DrugDiscovery.git
cd DrugDiscovery
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install -r requirements.txt
```

Dependencies are listed in [requirements.txt](../requirements.txt) at the repo root. The ML pipeline uses: `numpy`, `pandas`, `scipy`, `scikit-learn`, `lightgbm`, `joblib`, `rdkit`, `matplotlib`, `tqdm` (PyQt5/selfies in that file are for the GUI app, not the notebooks).

## Trained models

Best models `.joblib` files are not posted here due to their size. Generate them via notebooks `01 → 02 → 03`, or download from [OneDrive](https://1drv.ms/u/c/e8dfa526a80620bb/IQCuOso0IWgoSKAB_rlyJIV-AdUhLmeQ2q7e8UVmD9nhLsA?e=oEitfr) and extract into `data/models/`.

## References

1. **Gupta, A., et al.** — EGFR pIC50 prediction on ChEMBL data (2D descriptors, Extra Trees, random split). [PubMed](https://pubmed.ncbi.nlm.nih.gov/40718840/)
2. **Akinnusi, O. T., et al.** — EGFR activity prediction (ECFP + descriptors, Extra Trees, **scaffold split**). [Journal of Computer-Aided Molecular Design](https://link.springer.com/article/10.1007/s10822-026-00853-y)
3. **ChEMBL** — bioactive compound database. [https://www.ebi.ac.uk/chembl/](https://www.ebi.ac.uk/chembl/)
4. **RDKit** — cheminformatics and fingerprint generation. [https://www.rdkit.org/](https://www.rdkit.org/)
5. **Bemis, G. W.; Murcko, M. A.** — Murcko scaffold definition (used for splitting). *J. Med. Chem.* 1996. [https://pubmed.ncbi.nlm.nih.gov/8709122/](https://pubmed.ncbi.nlm.nih.gov/8709122/)

