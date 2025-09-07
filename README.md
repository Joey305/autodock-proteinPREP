# autodock-proteinPREP

A lightweight, interactive tool to prep protein receptors for **AutoDock Vina**.

* Reads **PDB** and **mmCIF**.
* Lists only **true HET groups** (ligands/ions/waters/sugars) by index for easy removal.
* Optional **chain** removal per file.
* **Collapses altLocs** to a single conformation per atom.
* Converts cleaned structures to **PDBQT** using whichever backend you have available.
* **CIF strict mode**: if any input is `.cif/.mmcif`, the final output folder contains **only** `.pdbqt` files.

> This repo **vendors** `AutoDockTools_py3/` (no submodules required). The environment installs it in editable mode automatically.

---

## What’s inside

```
3a_PDB2PDBQTbatch.py        # main interactive prep tool
AutoDockTools_py3/          # vendored ADT (Python 3 port)
Receptors/                  # example inputs (PDB + CIF mixed)
Receptors-Copy/             # tiny example set
environment.yml             # environment spec
```

Backends are auto-detected in this order:

1. **Meeko** (`mk_prepare_receptor.py`)
2. **AutoDockTools\_py3** (installed module from the vendored folder)
3. **AutoDockTools\_py3** (local checkout path with `PYTHONPATH`)
4. **MGLTools** legacy `prepare_receptor4.py` on PATH (optional fallback)

---

## Quick start

```bash
# 1) Clone
 git clone https://github.com/Joey305/autodock-proteinPREP.git
 cd autodock-proteinPREP

# 2) Create & activate the environment
 conda env create -f environment.yml
 conda activate vina

# 3) Run the prep tool
 python 3a_PDB2PDBQTbatch.py
```

> The environment installs `gemmi`, `meeko`, and `AutoDockTools_py3` (from the vendored folder) automatically.

---

## Example walkthroughs

### A) Small set: `Receptors-Copy/`

This folder includes a PDB and a CIF — so **CIF strict mode** will apply (only `.pdbqt` files will be kept in the output).

```text
📁 Available folders in current directory:
1. AutoDockTools_py3
2. Receptors
3. Receptors-Copy
🔍 Enter the number of the folder containing receptor structures (.pdb / .cif): 3

# Multiple files? You’ll be asked if you want batch mode.
🗂️ Multiple structures detected. Batch process with the SAME choices for all? [Y/n]: n   # choose per-file

── 5ND2.pdb ──
🔎 HET (ligands/ions/waters) in this file:
  1. HOH | count: ...
🔗 Chains in this file: C
🧽 Select HET indices to remove for THIS file (1,2,...) or 'all' or blank to skip: 1
🧩 Chains to remove for THIS file (comma-separated, blank to keep all):    # press Enter to keep C

── 7hg9.cif ──
🔎 HET ...
🔗 Chains in this file: A, B, C
🧽 Select HET indices ...: all
🧩 Chains to remove for THIS file: A,C

# Conversion runs; output goes to Receptors-Copy_PDBQT_Converted/
# Because there's a CIF, only .pdbqt files are kept in the output folder.
```

### B) Larger set: `Receptors/`

Mixed PDB + CIF; you can choose **batch mode** to apply the same HET/chain choices to all, or **per-file** to tailor each.

```text
🔍 Enter the number of the folder containing receptor structures (.pdb / .cif): 2
🗂️ Multiple structures detected. Batch process with the SAME choices for all? [Y/n]: y
🔎 HET ... (aggregated across all files)
🧽 Select HET indices to remove (e.g. 1,2,5) or 'all' to remove all shown: all
🧩 Remove specific chains for ALL files? (comma-separated, blank to skip): B,C
```

* Output folder: `Receptors_PDBQT_Converted/`
* Because this set contains a **CIF**, only `.pdbqt` files are kept in the output folder.

> Tip: Re-run in **per-file mode** when chain IDs differ among structures and you want different chain removal per file.

---

## Usage details

### HET selection

* The tool lists **true HET** residues only (ligands/ions/waters/sugars). Standard amino acids/nucleotides are never listed as HETs.
* Choose by **index** (e.g., `1,3,5`) or **`all`**. You can also add extra 3‑letter codes manually.

### Chain removal

* Provide chain IDs like `A,B` to drop whole chains. In per-file mode, you choose chains **per structure**.

### AltLoc handling

* AltLocs are **collapsed** automatically to one atom per name using highest occupancy (ties: `' '` > `A` > lexicographic). This prevents ADT’s alternate-location warnings.

### Outputs

* `<ChosenFolder>_PDBQT_Converted/` is created.
* For every input: `<basename>.converted.pdbqt`.
* If **only PDBs** were processed, the tool also writes:

  * `<basename>.clean.pdb` (cleaned PDB before conversion)
  * `HETATM_summary.txt` (HET counts across the set)
* If **any CIF** is present: output folder contains **only** `.pdbqt` (strict mode sweep).

---

## Installation notes

### Requirements

* Python 3.9–3.12 (3.11 recommended)
* Conda (or Python venv) and `pip`
* Linux/macOS/Windows (WSL recommended on Windows)

### What `environment.yml` does

* Installs `gemmi>=0.6.3` and `meeko>=0.5.0` via pip.
* Installs **vendored** `AutoDockTools_py3/` in editable mode (`-e ./AutoDockTools_py3`).

### Verify your setup

```bash
python -c "import gemmi; print('gemmi OK')"
python -c "import meeko; print('meeko OK')"          # optional but recommended
python -c "import AutoDockTools, MolKit; print('ADT_py3 OK')"
```

### Prefer pip‑installing ADT\_py3 instead of vendoring?

Remove the editable line from `environment.yml` and run:

```bash
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

The script will auto-detect it.

---

## Troubleshooting

**MolKit/AutoDockTools import error**

* You’re probably outside the env or the vendored install failed. Recreate the env, or install ADT\_py3 via pip from GitHub.

**`gemmi` not found**

* `pip install gemmi` (already handled by `environment.yml`).

**Thousands of altloc warnings in ADT**

* Safe to ignore; the script collapses altlocs first. If a particular structure still explodes with warnings, please open an issue and attach the file name.

**No receptor preparer found**

* Install at least one backend: `pip install meeko` **or** `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`.

**Windows tip**

* Use **WSL** for best results. Long paths and CRLF line endings can cause oddities in native shells.

---

## Tips

* Keep waters? Don’t select `HOH/WAT/H2O` when choosing HET indices.
* Different chain policies per structure? Choose **per-file** mode.
* Commit `HETATM_summary.txt` (PDB‑only runs) to document what was removed.

---

## License

MIT (recommended for tooling). If you use another license, add a `LICENSE` file.

---

## Acknowledgements

* **Gemmi** for fast, robust structure I/O.
* **Meeko** for modern AutoDock/Vina IO.
* **AutoDockTools\_py3** (Valdes‑Tresanco et al.) for the Python‑3 port of `prepare_receptor4.py`.
