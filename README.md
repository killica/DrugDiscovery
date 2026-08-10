# DrugDiscovery

Project for Computational intelligence course. Genetic algorithm for discovery of molecular structures with potential medical properties, using developed GUI application.

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

### Quick checklist


| Step                    | Command / action                                     |
| ----------------------- | ---------------------------------------------------- |
| Virtual env             | `python3 -m venv .venv && source .venv/bin/activate` |
| Install packages        | `pip install -r requirements.txt`                    |
| Start app               | `cd geneticAlgorithm && python main.py`              |
| pIC50 models (optional) | Run notebooks `data/01` → `02` → `03`                |

## 🖼️ Gallery

introafterLaunch

## 📸 Demo video

[https://github.com/user-attachments/assets/2181ebf9-822a-4bc1-841e-367bea711e4f](https://github.com/user-attachments/assets/2181ebf9-822a-4bc1-841e-367bea711e4f)