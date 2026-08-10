# DrugDiscovery

An interactive tool for **in silico** exploration of drug-like molecules. Instead of exhaustively testing large compound libraries in the lab, the application uses a **genetic algorithm** to evolve and rank structures against configurable fitness criteria — from drug-likeness (QED) to potency predictions from trained models.

Accompanying repository for a **Master's thesis at Faculty of Mathematics, University of Belgrade** on applying computational intelligence to computer-aided drug discovery.

**Thesis (PDF):** [thesis.pdf](thesis/build/thesis.pdf)

## 🔧 Setup

### 1. Prerequisites

- **Python 3.9 or newer** ([python.org](https://www.python.org/downloads/))
- `git` (to clone the repository)



### 2. Clone and enter the project

```bash
git clone <repository-url>
cd DrugDiscovery
```



### 3. Create a virtual environment (recommended)

**macOS / Linux:**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

**Windows (PowerShell):**

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```



### 4. Install dependencies

From the project root:

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

If `rdkit` fails to install via pip on your system, use conda instead:

```bash
conda create -n drugdiscovery python=3.11
conda activate drugdiscovery
conda install -c conda-forge rdkit pyqt matplotlib
pip install -r requirements.txt
```

(`conda` may already satisfy some packages; `pip` will skip duplicates or upgrade as needed.)

### 5. Run the GUI application

```bash
cd geneticAlgorithm
python main.py
```

The app works out of the box with the **QED** fitness mode. No trained models are required for that.

### 6. pIC50 fitness modes (optional)

Fitness modes **Random Forest**, **LightGBM**, and **Ridge** load pretrained models from `data/models/` (`*_best.joblib`). Those files are not in the repository — generate them once from the project root:

1. Open and run `data/01_preprocessing.ipynb`
2. Then `data/02_fingerprints.ipynb`
3. Then `data/03_models.ipynb`

After step 3, the GUI can use all four fitness modes.

Alternatively, these models are also available for download: (`.joblib` format): **[Download pretrained models (OneDrive)](https://1drv.ms/u/c/e8dfa526a80620bb/IQCuOso0IWgoSKAB_rlyJIV-AdUhLmeQ2q7e8UVmD9nhLsA?e=oEitfr)**

Extract the archive so that `data/models/` contains `ridge_best.joblib`, `lightgbm_best.joblib`, and `random_forest_best.joblib`.

### Quick checklist


| Step                    | Command / action                                                                         |
| ----------------------- | ---------------------------------------------------------------------------------------- |
| Virtual env             | `python3 -m venv .venv && source .venv/bin/activate`                                     |
| Install packages        | `pip install -r requirements.txt`                                                        |
| Start app               | `cd geneticAlgorithm && python main.py`                                                  |
| pIC50 models (optional) | - Run notebooks `data/01` → `02` → `03` OR Download from the OneDrive link - Restart app |




## 🖼️ Gallery

introafterLaunch

## 📸 Demo video

[https://github.com/user-attachments/assets/2181ebf9-822a-4bc1-841e-367bea711e4f](https://github.com/user-attachments/assets/2181ebf9-822a-4bc1-841e-367bea711e4f)