# autodock-proteinPREP

A lightweight, interactive tool to prep protein receptors for AutoDock Vina.
It reads **PDB** and **mmCIF**, lets you remove **true HET groups** (ligands/ions/waters/sugars) and/or **specific chains**, collapses **altLocs**, and converts everything to **PDBQT** using whichever backend you have available (Meeko or AutoDockTools).

---

## Why this exists

* I didn’t want to memorize half a dozen prep commands for every receptor source.
* I needed a single script that works in different environments (with or without Meeko / ADT).
* I wanted fast, repeatable cleanup: pick HETs by index, optionally drop chains, get PDBQT out.

---

## Quick start (with submodules)

```bash
# 1) Clone WITH SUBMODULES (needed for AutoDockTools_py3)
git clone --recurse-submodules https://github.com/Joey305/autodock-proteinPREP.git
cd autodock-proteinPREP

# 2) Create & activate the environment from this repo
conda env create -f environment.yml
conda activate vina

# 3) Run the receptor prep tool
python 3a_PDB2PDBQTbatch.py
```

> Already cloned without submodules?
>
> ```bash
> git submodule update --init --recursive
> ```

---

## What it does

* Accepts **.pdb / .ent / .cif / .mmcif** from a folder you choose.
* Scans and lists only **true HET** residue names (ligands/ions/waters/sugars).
* Lets you select HETs by **index** (e.g., `1,3,5`) or **`all`**.
* **Optional chain removal** by chain ID (e.g., `A,B`).
* **AltLoc collapse**: keeps one atom per name, preferring highest occupancy (ties: `' '` > `A` > lexicographic).
* Converts cleaned structures to **PDBQT** using the first available backend:

  1. **Meeko** `mk_prepare_receptor.py`
  2. **AutoDockTools\_py3** (module installed)
  3. **AutoDockTools\_py3** (local submodule checkout)
  4. **MGLTools** `prepare_receptor4.py` on PATH (legacy)
* **CIF strict mode**: if **any** input is CIF/mmCIF, the output folder contains **only** `.pdbqt` files (no `.clean.pdb`, no summary).

---

## Installation

**Requirements**

* Python **3.9–3.12** (3.11 recommended)
* `conda` or `pip`
* Linux/macOS or Windows (WSL recommended)

### Environment (from this repo)

The included `environment.yml` installs:

* `gemmi` (structure IO)
* `meeko` (modern AutoDock IO)
* `AutoDockTools_py3` from the **submodule** (editable install `-e ./AutoDockTools_py3`)

Create it like this:

```bash
conda env create -f environment.yml
conda activate vina
```

### Alternative: no submodules

If you don’t want to use submodules, you can remove the editable line in `environment.yml` and install ADT\_py3 directly from GitHub:

```bash
pip install gemmi meeko
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

The script auto-detects which backend you have and uses that.

---

## Usage

Put your receptor files under a folder like `Receptors/`, then run:

```bash
python 3a_PDB2PDBQTbatch.py
```

You’ll see:

1. A prompt to choose a folder containing structures.
2. If multiple files are found, you can **batch** process with the same choices, or handle them **per-file** with different HET/chain selections.
3. A numbered list of **HET** groups (true ligands/ions/waters/sugars) to remove.
4. An optional prompt to remove **chains** by ID.

### Modes

* **Batch mode**: one HET/chain selection is applied to **all** structures.
* **Per-file mode**: you make HET/chain choices **independently** per structure.

### AltLoc handling

Alternate locations are collapsed automatically to one conformation per atom using highest occupancy (tie-breakers: `' '` > `A` > lexicographic). This prevents the “thousands of alternate location atoms” warnings in ADT.

### CIF strict mode

If **any** input file is `.cif/.mmcif`, the output folder will contain **only** `.pdbqt`.
No `.clean.pdb` and no summary file are kept in that case.

### Outputs

* Output folder: `<ChosenFolder>_PDBQT_Converted/`
* For every input: `<basename>.converted.pdbqt`
* If **only PDBs** were processed, the tool also writes:

  * `<basename>.clean.pdb` (cleaned input)
  * `HETATM_summary.txt` (HET counts across the set)

---

## Troubleshooting

**“MolKit / AutoDockTools import error”**
Make sure you either:

* cloned with submodules and created the env from `environment.yml`, or
* installed `AutoDockTools_py3` via `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`.

**“gemmi not found”**
Install it: `pip install gemmi` (already in `environment.yml`).

**ADT warnings about altlocs**
These are safe to ignore; the script already collapses altlocs before conversion. If a particular PDB/mmCIF still triggers huge warnings, open an issue with the file.

**No receptor preparer found**
Install at least one backend:

* `pip install meeko` (recommended), or
* `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`, or
* ensure legacy MGLTools is on your PATH.

---

## Submodule management (for contributors)

Update ADT\_py3 to its latest upstream:

```bash
git submodule update --remote --merge
git add AutoDockTools_py3
git commit -m "Update AutoDockTools_py3 submodule"
```

If you prefer SSH for submodules:

```bash
git config -f .gitmodules submodule.AutoDockTools_py3.url git@github.com:Valdes-Tresanco-MS/AutoDockTools_py3.git
git submodule sync --recursive
git submodule update --init --recursive
```

---

## Tips

* Want to keep waters? Don’t select `HOH/WAT/H2O` when choosing HET indices.
* Chain IDs can differ per structure. Use **per-file** mode to tailor removals.
* Commit the `HETATM_summary.txt` (when PDB-only) to record what was removed.

---

## License

MIT (recommended for tooling). If you adopt another license, drop it in `LICENSE`.

---

## Acknowledgements

* **Gemmi** for reliable structure IO.
* **Meeko** for modern AutoDock/Vina IO.
* **AutoDockTools\_py3** (Valdes-Tresanco et al.) for the Python-3 port of `prepare_receptor4.py`.
