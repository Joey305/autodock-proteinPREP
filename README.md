# autodock-proteinPREP

A lightweight, interactive **and** headless tool to prep protein receptors for **AutoDock Vina**.

* Reads **PDB** and **mmCIF**
* Lists only **true HET groups** (ligands/ions/waters/sugars) by index
* Optional **chain** removal (batch or per‑file)
* **Collapses altLocs** to a single conformation per atom
* Converts to **PDBQT** via Meeko / AutoDockTools (ADT) / legacy MGLTools
* **PDBQT‑only output** by default (final folder contains only `.pdbqt`)

This repository **vendors** `AutoDockTools_py3/` (no submodule steps). The environment installs it in editable mode automatically.

---

## Quick start

```bash
# 1) Clone
 git clone https://github.com/Joey305/autodock-proteinPREP.git
 cd autodock-proteinPREP

# 2) Create & activate the environment
 conda env create -f environment.yml
 conda activate vina

# 3) Run interactively (heads‑up display)
 python 3a_PDB2PDBQTbatch.py
```

You’ll pick a folder, review HETs and chains, choose removals, and the tool writes PDBQT files to `<Folder>_PDBQT_Converted/`.

---

## What it does

* Accepts **.pdb / .ent / .cif / .mmcif** from a chosen folder
* Detects **true HET** residues only (never lists standard amino acids or nucleotides)
* Lets you remove HETs by **index** (or `all`) and optionally remove **chains** by ID
* **AltLoc** policy: collapses to one atom per name using highest occupancy (tie: `' '` > `A` > lexicographic)
* **Backends (auto‑detect, in order):**

  1. **Meeko** (`mk_prepare_receptor.py`)
  2. **AutoDockTools\_py3** (installed module)
  3. **AutoDockTools\_py3** (local vendored path)
  4. **MGLTools** legacy `prepare_receptor4.py` on PATH
* **PDBQT‑only**: by default the final output folder contains only `.pdbqt` files (non‑PDBQT artifacts are swept)

---

## Installation

**Requirements**

* Python **3.9–3.12** (3.11 recommended)
* Conda (or another venv) and `pip`
* Linux/macOS or Windows (WSL recommended on Windows)

**Environment**

The included `environment.yml` installs:

* `gemmi` (structure I/O)
* `meeko` (modern AutoDock/Vina I/O)
* vendored **AutoDockTools\_py3** in editable mode (`-e ./AutoDockTools_py3`)

Create the environment:

```bash
conda env create -f environment.yml
conda activate vina
```

Verify:

```bash
python -c "import gemmi; print('gemmi OK')"
python -c "import meeko; print('meeko OK')"            # optional but recommended
python -c "import AutoDockTools, MolKit; print('ADT_py3 OK')"
```

> Prefer to install ADT\_py3 from GitHub instead of vendoring? Replace the editable line in `environment.yml` with:
> `git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`

---

## Example walkthroughs (from the repo)

### A) Small set — `Receptors-Copy/`

This folder includes a PDB and a CIF. Run:

```bash
python 3a_PDB2PDBQTbatch.py
```

* Select **`Receptors-Copy`** at the prompt
* Choose **per‑file** mode if you want different chain/HET policies per file
* Because this set contains a **CIF**, the output folder will end with **only `.pdbqt`** (no summary / no `.clean.pdb`)

### B) Larger set — `Receptors/`

Run interactively and choose **batch** to apply the same choices across all files, or **per‑file** to tailor each.

> Tip: If chain IDs differ and you need per‑structure choices, pick **per‑file** mode.

---

## Headless / automation

Everything you can do interactively can be done headlessly (no prompts). Defaults in headless mode: **batch**, **remove all HETs**, **keep all chains**, **PDBQT‑only on**.

**Batch, remove all HETs, keep all chains (defaults):**

```bash
python 3a_PDB2PDBQTbatch.py \
  --folder Receptors \
  --headless
```

**Batch, remove specific HET codes and chains:**

```bash
python 3a_PDB2PDBQTbatch.py \
  --folder Receptors \
  --headless --mode batch \
  --remove-het HOH,EDO,SO4 \
  --remove-chains B,C
```

**Batch, pick HETs by aggregated indices (from the scan shown in logs):**

```bash
python 3a_PDB2PDBQTbatch.py \
  --folder Receptors \
  --headless --mode batch \
  --remove-het-indices 1,3,5
```

**Per‑file headless via JSON config:**

`config.json`

```json
{
  "7hg9.cif": { "remove_het": "all", "remove_chains": ["A", "C"] },
  "8BB2_cleaned.pdb": { "remove_het": ["HOH", "EDO"], "remove_chains": [] }
}
```

Run:

```bash
python 3a_PDB2PDBQTbatch.py \
  --folder Receptors \
  --headless --mode per-file \
  --per-file-config config.json
```

**Limit to a subset of files:**

```bash
python 3a_PDB2PDBQTbatch.py \
  --folder Receptors \
  --files 7hg9.cif 8BB2_cleaned.pdb \
  --headless
```

**Pin the backend, adjust altLoc policy, and allow artifacts:**

```bash
python 3a_PDB2PDBQTbatch.py \
  --folder Receptors \
  --headless \
  --backend meeko \
  --altloc collapse \
  --no-pdbqt-only \
  --write-summary \
  --keep-clean-pdb
```

### CLI reference (most useful)

```
--folder FOLDER                Input folder containing .pdb/.cif
--files F1 [F2 ...]            Optional subset of files to process (relative to --folder)

--headless                     Run without prompts
--mode {batch,per-file}        Headless selection mode (default: batch)

--remove-het ALL|CODES         Comma list (e.g. HOH,EDO,SO4) or 'all' (default in headless: all)
--remove-het-indices IDS       Comma list of indices from aggregated HET scan (batch only)
--remove-chains IDS            Comma list of chain IDs (e.g. A,B) (default headless: none)
--per-file-config JSON         JSON mapping filename -> {remove_het, remove_chains}

--backend {auto,meeko,adt,mgl} Backend preference (default: auto)
--altloc {collapse,all}        AltLoc handling (default: collapse)

--output-dir PATH              Output directory (default: <folder>_PDBQT_Converted)
--pdbqt-only / --no-pdbqt-only Keep only .pdbqt in final folder (default: on)
--write-summary                Write HET summary (suppressed if --pdbqt-only)
--keep-clean-pdb               Save cleaned PDBs (suppressed if --pdbqt-only)
```

**Headless defaults**: `--mode batch`, `--remove-het all`, **keep all chains**, `--pdbqt-only`.

---

## Outputs

* Output folder: `<Folder>_PDBQT_Converted/`
* For every input structure: `<basename>.converted.pdbqt`
* If you *disable* PDBQT‑only (with `--no-pdbqt-only`):

  * `.clean.pdb` copies can be kept (`--keep-clean-pdb`)
  * `HETATM_summary.txt` can be written (`--write-summary`; PDB inputs only)

---

## Troubleshooting

**MolKit / AutoDockTools import error**

* Ensure you’re inside the `vina` env created from `environment.yml` **or** install ADT\_py3 from GitHub:
  `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`

**`gemmi` not found**

* `pip install gemmi` (already handled by the environment file)

**Huge ADT altloc warnings**

* Safe to ignore; the script collapses altlocs before conversion. If a specific file still floods warnings, open an issue with the filename.

**No receptor preparer found**

* Install at least one backend:

  * `pip install meeko` (recommended), or
  * `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`, or
  * ensure MGLTools’ `prepare_receptor4.py` is on PATH.

**Windows tips**

* Prefer **WSL**. Native Windows shells (path length, CRLF) can be finicky for structural pipelines.

---

## Tips

* Keep waters? Don’t select `HOH/WAT/H2O` when choosing HET indices.
* Different chain policies per structure? Use **per‑file** mode (interactive or headless).
* For audit trails, disable PDBQT‑only and enable `--write-summary` and `--keep-clean-pdb`.

---

## License

MIT (recommended for tooling). If you use another license, add a `LICENSE` file.

---

## Acknowledgements

* **Gemmi** for fast, robust structure I/O
* **Meeko** for modern AutoDock/Vina I/O
* **AutoDockTools\_py3** (Valdes‑Tresanco et al.) for the Python‑3 port of `prepare_receptor4.py`
