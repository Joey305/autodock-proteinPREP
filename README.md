# Vina‑ProteinPrep

Batch and per‑file receptor preparation for AutoDock Vina. This repo provides a single interactive script that:

* Reads **PDB** and **mmCIF** files.
* Detects and lists only **true HET groups** (ligands/ions/waters/sugars) with indexed selection.
* Optionally removes **specific chains**.
* **Collapses alternate locations** (altLoc) to one conformation per atom.
* Converts cleaned structures to **PDBQT** using whichever tool you have available, in this order:

  1. **Meeko** (`mk_prepare_receptor.py`)
  2. **AutoDockTools\_py3** (installed module)
  3. **AutoDockTools\_py3** (local checkout in the repo)
  4. **Legacy MGLTools** `prepare_receptor4.py` on PATH
* **CIF strict mode:** if any input is `.cif/.mmcif`, the output folder will contain **only** `.pdbqt` files (no extra artifacts).

> The goal is a friction‑free, shareable workflow that adapts to different environments.

---

## Table of contents

* [Quick start](#quick-start)
* [Installation](#installation)

  * [Option A — Meeko only (simplest)](#option-a--meeko-only-simplest)
  * [Option B — AutoDockTools\_py3 (ADT, Python 3)](#option-b--autodocktools_py3-adt-python-3)
  * [Option C — Legacy MGLTools (Python 2, optional fallback)](#option-c--legacy-mgltools-python-2-optional-fallback)
  * [Environment files](#environment-files)
  * [Verify your setup](#verify-your-setup)
* [Usage](#usage)

  * [Modes: batch vs per‑file](#modes-batch-vs-per-file)
  * [HET selection](#het-selection)
  * [Chain removal](#chain-removal)
  * [AltLoc handling](#altloc-handling)
  * [CIF strict mode](#cif-strict-mode)
  * [Outputs](#outputs)
* [How it works](#how-it-works)
* [Troubleshooting](#troubleshooting)
* [Tips](#tips)
* [Contributing](#contributing)
* [Acknowledgements](#acknowledgements)
* [License](#license)

---

## Quick start

```bash
# 1) Clone
git clone https://github.com/<you>/Vina-ProteinPrep.git
cd Vina-ProteinPrep

# 2) Create & activate the environment from this repo
conda env create -f environment.yml
conda activate vina

# 3) (Optional) If you didn't vendor AutoDockTools_py3/ in this repo,
#    install the ADT_py3 backend from GitHub instead:
# pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3

# 4) Run the receptor prep tool
python 3a_PDB2PDBQTbatch.py
```

---

## Installation

**Requirements**

* Python **3.9–3.12** (3.11 recommended)
* `pip` or `conda`
* OS: Linux, macOS, or Windows (PowerShell/WSL recommended on Windows)

You can use **Meeko**, **AutoDockTools\_py3**, **or** legacy **MGLTools**. The script auto‑detects and uses the first available.

### Option A — Meeko only (simplest)

Meeko is the actively maintained I/O stack for AutoDock/Vina.

```bash
pip install meeko gemmi
```

No further setup required. The script will use `mk_prepare_receptor.py`.

### Option B — AutoDockTools\_py3 (ADT, Python 3)

Install the maintained Python‑3 port of ADT to get `prepare_receptor4.py` without Python 2:

```bash
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

**Alternative (vendored checkout):** If this repo contains an `AutoDockTools_py3/` folder, you can install it locally:

```bash
pip install -e ./AutoDockTools_py3
```

> Editable install is convenient for development; for production, prefer `pip install` from the GitHub URL above.

### Option C — Legacy MGLTools (Python 2, optional fallback)

If you must use classic MGLTools, put `prepare_receptor4.py` on your PATH (e.g., via a separate Python 2 env) and the script will fall back to it automatically. This is optional.

### Environment files

A minimal `environment.yml` you can ship with this repo:

```yaml
ame: vina
channels: [conda-forge, defaults]
dependencies:
  - python=3.11
  - pip
  - pip:
      - gemmi>=0.6.3
      - meeko>=0.5.0
```

Usage:

```bash
conda env create -f environment.yml
conda activate vina
```

Then optionally install ADT\_py3:

```bash
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

### Verify your setup

```bash
# Gemmi present?
python -c "import gemmi; print('gemmi OK')"

# Meeko present (optional)?
python -c "import meeko; print('meeko OK')"  # prints if installed

# ADT_py3 present (optional)?
python -c "import AutoDockTools, MolKit, mglutil; print('ADT_py3 OK')"
```

---

## Usage

Place your receptor structures under a folder like `Receptors/` (PDB and/or mmCIF).

Start the interactive prep:

```bash
python 3a_PDB2PDBQTbatch.py
```

You’ll be prompted to pick a folder. If multiple structures are found, you’ll be asked whether to **batch** process (same choices for all) or go **per‑file**.

### Modes: batch vs per‑file

* **Batch mode**: one set of HET removals and chain removals is applied to all structures.
* **Per‑file mode**: you’ll review **each file**, select HET indices and chains to remove independently.

### HET selection

The script scans for **true HET groups** only:

* Uses Gemmi entity typing (NonPolymer/Water) and excludes standard amino acids and nucleotides.
* Treats common glycans (e.g., `NAG`, `BMA`, `MAN`, `GAL`, `FUC`, `NDG`) as HETs.
* Presents a numbered list (e.g., `1: HOH`, `2: SO4`, `3: EDO`, …). Choose indices (e.g., `1,3,5`) or `all`.
* You can add extra 3‑letter codes manually.

### Chain removal

Optionally remove chains by ID (e.g., `A,B`). In **per‑file mode**, chain choices are independent per structure.

### AltLoc handling

Alternate locations are **collapsed automatically** to a single conformation per atom using the policy:

* highest **occupancy** wins; ties break by altLoc `' '` then `'A'`, then lexicographic.
  This prevents ADT warnings about thousands of alternate‑location atoms.

### CIF strict mode

If **any** input file is `.cif/.mmcif`, the output folder will contain **only** `.pdbqt` files. No `.clean.pdb` or summary files are kept.

### Outputs

* Output folder: `<ChosenFolder>_PDBQT_Converted/`
* For each input: `<basename>.converted.pdbqt`
* When **only PDBs** are processed, the script also writes:

  * `<basename>.clean.pdb` (post‑cleanup, pre‑conversion)
  * `HETATM_summary.txt` (counts of detected HETs across the set)

---

## How it works

1. **Scan** structures with Gemmi, listing HET residue names and available chains.
2. **Select** HETs by index (`all` supported) and optional extra residue codes.
3. **Optionally remove chains** by ID.
4. **Clean** each structure: drop selected HET residues and chains, **collapse altLocs**, and write a clean PDB via a version‑agnostic writer (works regardless of Gemmi version).
5. **Convert to PDBQT** using the first available preparer:

   * Meeko → `mk_prepare_receptor.py`
   * ADT\_py3 (installed) → `python -m AutoDockTools.Utilities24.prepare_receptor4`
   * ADT\_py3 (local checkout) → direct script with `PYTHONPATH`
   * MGLTools (legacy) → `prepare_receptor4.py`
6. **CIF strict sweep**: if any CIF was processed, remove any non‑`.pdbqt` artifacts from the output folder.

---

## Troubleshooting

**“MolKit not found” or ADT import errors**

* Install ADT\_py3: `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`.
* Or ensure the repo’s `AutoDockTools_py3/` is installed: `pip install -e ./AutoDockTools_py3`.

**Gemmi not installed**

* `pip install gemmi`

**Thousands of altLoc warnings from ADT**

* The script collapses altLocs automatically. If you still see warnings, your source file may contain unusual altLoc labels; open an issue with the problematic PDB/mmCIF.

**No receptor preparer found**

* Install one of: Meeko, ADT\_py3, or add MGLTools to PATH. The script prints what it tried.

**Windows notes**

* Prefer **WSL** or PowerShell. Long paths and permissions may differ; ensure the repo path has write access.

---

## Tips

* Want to keep waters but drop ligands? Select indices that exclude `HOH`/`WAT`.
* Use **per‑file mode** when chain IDs differ between structures.
* Keep the `HETATM_summary.txt` committed (when PDB‑only) to document data cleaning choices.

---

## Contributing

Issues and PRs are welcome! If you add features (e.g., CLI flags for non‑interactive runs), please include tests or clear examples.

---

## Acknowledgements

* **Gemmi** for fast, robust structural I/O.
* **Meeko** for modern AutoDock/Vina I/O.
* **AutoDockTools\_py3** for the Python‑3 port of `prepare_receptor4.py`.

---

## License

Choose a license (e.g., MIT) and add a `LICENSE` file. If you’re unsure, MIT is a good default for tooling like this.
