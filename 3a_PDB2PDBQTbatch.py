#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3a_PDB2PDBQTbatch.py
Interactive receptor prep for AutoDock Vina with:
- PDB + mmCIF input support
- True-HET detection (ligands/ions/waters/sugars)
- Optional chain removal (batch or per-file)
- AltLoc collapse to a single conformation per atom
- Conversion to PDBQT via Meeko / ADT_py3 / legacy MGLTools
- FINAL OUTPUT: **only .pdbqt** files (always)

Requires: gemmi  (pip install gemmi)
Optional: meeko  (pip install meeko)
Optional: ADT_py3 (pip install -e ./AutoDockTools_py3 OR pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3)
"""

import os
import sys
import shutil
import subprocess
from collections import Counter, defaultdict
import importlib.util

# =============== Global policy ===============
# Always keep ONLY .pdbqt files in the final output directory
ALWAYS_PDBQT_ONLY = True

# =============== Dependencies ===============
try:
    import gemmi
except ImportError:
    print("❌ Requires 'gemmi'. Install with: pip install gemmi")
    sys.exit(1)

# =============== Helpers: module presence & backends ===============
def has_module(modname: str) -> bool:
    return importlib.util.find_spec(modname) is not None

def resolve_receptor_preparer():
    """
    Returns (kind, base_cmd_list, extra_env_or_None)
    kind ∈ {'meeko','adt_module','adt_local','mgltools',None}
    Priority:
      1) Meeko (mk_prepare_receptor.py)
      2) ADT_py3 installed as module
      3) ADT_py3 local checkout in ./AutoDockTools_py3 (PYTHONPATH)
      4) MGLTools legacy (prepare_receptor4.py on PATH)
    """
    meeko = shutil.which("mk_prepare_receptor.py")
    if meeko:
        return "meeko", [meeko], None

    if has_module("AutoDockTools.Utilities24.prepare_receptor4"):
        return "adt_module", [sys.executable, "-m", "AutoDockTools.Utilities24.prepare_receptor4"], None

    local_script = os.path.join("AutoDockTools_py3", "AutoDockTools", "Utilities24", "prepare_receptor4.py")
    if os.path.exists(local_script):
        env = os.environ.copy()
        adt_root = os.path.abspath("AutoDockTools_py3")
        env["PYTHONPATH"] = adt_root + (os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else "")
        return "adt_local", [sys.executable, local_script], env

    mgl = shutil.which("prepare_receptor4.py")
    if mgl:
        return "mgltools", [mgl], None

    return None, None, None

# =============== Constants for HET classification ===============
STD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "SEC","PYL","MSE"
}
STD_NT = {"A","C","G","T","U","I","DA","DC","DG","DT","DI","RA","RC","RG","RU"}
WATER_NAMES = {"HOH","WAT","H2O"}
COMMON_SUGARS = {"NAG","BMA","MAN","GAL","FUC","NDG"}  # treat as HET

# =============== HET classification ===============
def is_true_het(res: gemmi.Residue, chain=None) -> bool:
    name = res.name.strip().upper()
    if name in WATER_NAMES:
        return True
    et = res.entity_type
    if et in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Water):
        return True
    if name in STD_AA or name in STD_NT:
        return False
    if name in COMMON_SUGARS:
        return True
    if res.het_flag != ' ':
        return True
    return False

# =============== AltLoc policy ===============
def iter_atoms_with_altloc_policy(res, policy="collapse"):
    """
    'collapse': choose a single atom per atom name using highest occupancy;
                ties broken by altLoc ' ' then 'A' then lexicographic.
    'all': yield all atoms (no filtering).
    """
    if policy == "all":
        for a in res:
            yield a
        return

    groups = {}
    for a in res:
        key = a.name.strip()
        alt = getattr(a, "altloc", "") or ""
        occ = getattr(a, "occ", 1.0)
        prev = groups.get(key)
        if prev is None:
            groups[key] = (a, alt, occ)
        else:
            _, p_alt, p_occ = prev
            better = (occ > p_occ) or (occ == p_occ and (alt or " ") < (p_alt or " "))
            if better:
                groups[key] = (a, alt, occ)
    for a, _, _ in groups.values():
        yield a

# =============== PDB writing (version-agnostic) ===============
def format_atom_name(atom_name: str) -> str:
    n = atom_name.strip()
    if len(n) >= 4:
        return n[:4]
    return f"{n:>4}"

def derive_element(atom) -> str:
    try:
        elem = atom.element.name.strip()
        if elem and elem != "X":
            return elem
    except Exception:
        pass
    letters = "".join([c for c in atom.name if c.isalpha()])
    if not letters:
        return ""
    if len(letters) >= 2 and letters[1].islower():
        return letters[:2].title()
    return letters[0].upper()

def pdb_write_manual(structure: gemmi.Structure, out_pdb_path: str):
    serial = 1
    lines = []
    wrote_any = False

    for mi, model in enumerate(structure):
        for chain in model:
            chain_id = (chain.name[:1] if chain.name else "A")
            chain_wrote = False
            last_resname, last_resseq, last_icode = "UNK", 0, " "

            for res in chain:
                resname = res.name.strip().upper()
                try:
                    resseq = res.seqid.num
                    icode = res.seqid.icode if res.seqid.icode and res.seqid.icode != "\x00" else " "
                except Exception:
                    resseq, icode = 0, " "
                record = "HETATM" if is_true_het(res, chain) else "ATOM  "
                for atom in iter_atoms_with_altloc_policy(res, policy="collapse"):
                    x, y, z = atom.pos.x, atom.pos.y, atom.pos.z
                    occ = getattr(atom, "occ", 1.00)
                    bfac = getattr(atom, "b_iso", 0.00)
                    altloc = atom.altloc if getattr(atom, "altloc", "") else " "
                    atname = format_atom_name(atom.name)
                    element = derive_element(atom)
                    line = (
                        f"{record:<6}{serial:>5} "
                        f"{atname:<4}{altloc:1}"
                        f"{resname:>3} "
                        f"{chain_id:1}{resseq:>4}{icode:1}"
                        f"   "
                        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                        f"{occ:>6.2f}{bfac:>6.2f}"
                        f"          {element:>2}  "
                    )
                    lines.append(line)
                    serial += 1
                    wrote_any = True
                    chain_wrote = True
                    last_resname, last_resseq, last_icode = resname, resseq, icode

            if chain_wrote:
                ter = f"{'TER':<6}{serial:>5}      {last_resname:>3} {chain_id:1}{last_resseq:>4}{last_icode:1}"
                lines.append(ter)
                serial += 1

    lines.append("END")
    with open(out_pdb_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if not wrote_any:
        raise RuntimeError("No atoms written — structure empty after filtering.")

def write_pdb_compat(structure: gemmi.Structure, out_pdb_path: str):
    if hasattr(gemmi, "write_minimal_pdb"):
        try:
            gemmi.write_minimal_pdb(structure, out_pdb_path)
            return
        except Exception:
            pass
    pdb_write_manual(structure, out_pdb_path)

# =============== IO helpers ===============
def read_structure(path: str) -> gemmi.Structure:
    return gemmi.read_structure(path)

def is_cif_name(filename: str) -> bool:
    ext = os.path.splitext(filename.lower())[1]
    return ext in {".cif", ".mmcif"}

# =============== Scanning helpers ===============
def scan_file_for_hets_and_chains(path: str):
    """Return (Counter{name->count}, set(chains)) for a single file."""
    hets = Counter()
    chains = set()
    st = read_structure(path)
    for model in st:
        for chain in model:
            chains.add(chain.name)
            for res in chain:
                if is_true_het(res, chain):
                    hets[res.name.strip().upper()] += 1
    return hets, chains

def scan_all(files, folder):
    """Aggregate HETs across all files + chains per file."""
    het_counter = Counter()
    het_presence = defaultdict(set)
    chains_per_file = defaultdict(set)
    for fname in files:
        fpath = os.path.join(folder, fname)
        try:
            h, chs = scan_file_for_hets_and_chains(fpath)
            for k, v in h.items():
                het_counter[k] += v
                het_presence[k].add(fname)
            chains_per_file[fname] = chs
        except Exception as e:
            print(f"⚠️ Skipping {fname}: {e}")
    return het_counter, het_presence, chains_per_file

# =============== Cleaning step ===============
def clean_structure_to_pdb(in_path, out_pdb_path, residues_to_drop, chains_to_drop):
    st = gemmi.read_structure(in_path)
    new_st = gemmi.Structure()
    new_st.cell = st.cell
    new_st.spacegroup_hm = st.spacegroup_hm
    new_st.name = st.name

    res_drop = {r.upper() for r in residues_to_drop}
    chain_drop = set(chains_to_drop)

    for mi, model in enumerate(st):
        new_model = gemmi.Model(str(mi + 1))
        for chain in model:
            if chain.name in chain_drop:
                continue
            new_chain = gemmi.Chain(chain.name)
            for res in chain:
                if is_true_het(res, chain) and res.name.strip().upper() in res_drop:
                    continue
                new_chain.add_residue(res)
            if len(new_chain) > 0:
                new_model.add_chain(new_chain)
        if len(new_model) > 0:
            new_st.add_model(new_model)

    write_pdb_compat(new_st, out_pdb_path)

# =============== Main ===============
def main():
    # list folders
    print("\n📁 Available folders in current directory:")
    folders = [f for f in os.listdir() if os.path.isdir(f)]
    for i, folder in enumerate(folders, 1):
        print(f"{i}. {folder}")
    choice = input("🔍 Enter the number of the folder containing receptor structures (.pdb / .cif): ")
    try:
        folder_index = int(choice) - 1
        selected_folder = folders[folder_index]
    except (ValueError, IndexError):
        raise ValueError("❌ Invalid folder selection.")

    # files
    valid_exts = {".pdb", ".ent", ".cif", ".mmcif"}
    files = sorted(
        f for f in os.listdir(selected_folder)
        if os.path.isfile(os.path.join(selected_folder, f))
        and os.path.splitext(f.lower())[1] in valid_exts
    )
    if not files:
        print("❌ No .pdb / .cif / .mmcif files found in the selected folder.")
        sys.exit(1)

    # Ask batch vs per-file if multiple
    batch_mode = True
    if len(files) > 1:
        ans = input("\n🗂️ Multiple structures detected. Batch process with the SAME choices for all? [Y/n]: ").strip().lower()
        batch_mode = (ans in ("", "y", "yes"))

    # Resolve preparer
    prep_kind, prep_base, prep_env = resolve_receptor_preparer()
    if not prep_kind:
        raise RuntimeError(
            "No receptor preparer found. Install Meeko (`pip install meeko`) or ADT_py3 "
            "(`pip install -e ./AutoDockTools_py3` or `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`) "
            "or ensure MGLTools' prepare_receptor4.py is on PATH."
        )

    # Output / temp
    output_dir = selected_folder + "_PDBQT_Converted"
    os.makedirs(output_dir, exist_ok=True)
    temp_cleaned_dir = ".temp_cleaned_pdbs"
    os.makedirs(temp_cleaned_dir, exist_ok=True)

    # Info banner
    if ALWAYS_PDBQT_ONLY:
        print("\n🧾 PDBQT-only mode: final output folder will contain ONLY .pdbqt files.\n")

    # Gather choices
    per_file_choices = {}

    if batch_mode:
        print(f"\n📦 Found {len(files)} structure file(s). Scanning for HET ligands/waters and chains...")
        het_counter, het_presence, chains_per_file = scan_all(files, selected_folder)

        het_names_sorted = sorted(het_counter.items(), key=lambda x: (-x[1], x[0]))
        print("\n🔎 HET (ligands/ions/waters) detected across all files:")
        if not het_names_sorted:
            print("  (none)")
        else:
            for idx, (name, count) in enumerate(het_names_sorted, 1):
                files_here = sorted(het_presence[name])
                print(f"  {idx:>3}. {name:<6} | count: {count:<4} | files: {', '.join(files_here)}")

        print("\n🔗 Chains detected (per file):")
        for fname in files:
            chs = sorted(chains_per_file.get(fname, []))
            print(f"  • {fname}: {', '.join(chs) if chs else '(none)'}")

        residues_to_remove = set()
        if het_names_sorted:
            sel = input("\n🧽 Select HET indices to remove (e.g. 1,2,5) or 'all' to remove all shown. Leave blank to skip: ").strip().lower()
            if sel == "all":
                residues_to_remove = {name for name, _ in het_names_sorted}
            elif sel:
                try:
                    chosen = {int(x.strip()) for x in sel.split(",") if x.strip()}
                    for i in chosen:
                        if 1 <= i <= len(het_names_sorted):
                            residues_to_remove.add(het_names_sorted[i-1][0])
                        else:
                            print(f"  ⚠️ Ignoring out-of-range index: {i}")
                except ValueError:
                    print("  ⚠️ Invalid index list, skipping automatic HET removal.")

        extra_res = input("➕ Enter any EXTRA 3-letter residue names to remove (comma-separated, optional): ").strip()
        if extra_res:
            residues_to_remove |= {r.strip().upper() for r in extra_res.split(",") if r.strip()}
        if residues_to_remove:
            print("🧼 Will remove HET residues (for all files):", ", ".join(sorted(residues_to_remove)))
        else:
            print("ℹ️ No HET residue names selected for removal.")

        chain_resp = input("🧩 Remove specific chains for ALL files? Enter chain IDs (comma-separated, leave blank to skip): ").strip()
        chains_to_remove_global = {c.strip() for c in chain_resp.split(",") if c.strip()} if chain_resp else set()
        if chains_to_remove_global:
            print("✂️ Will remove chains (for all files):", ", ".join(sorted(chains_to_remove_global)))
        else:
            print("ℹ️ No chains selected for removal.")

        for fname in files:
            per_file_choices[fname] = (set(residues_to_remove), set(chains_to_remove_global))

    else:
        print(f"\n📦 Found {len(files)} structure file(s). Per-file selection enabled.\n")
        for fname in files:
            fpath = os.path.join(selected_folder, fname)
            print(f"── {fname} ─────────────────────────────────────────────")
            try:
                hets, chains = scan_file_for_hets_and_chains(fpath)
            except Exception as e:
                print(f"⚠️ Skipping (scan failed): {e}")
                continue

            het_sorted = sorted(hets.items(), key=lambda x: (-x[1], x[0]))
            print("🔎 HET (ligands/ions/waters) in this file:")
            if not het_sorted:
                print("  (none)")
            else:
                for idx, (name, count) in enumerate(het_sorted, 1):
                    print(f"  {idx:>3}. {name:<6} | count: {count:<4}")

            print("🔗 Chains in this file:", ", ".join(sorted(chains)) if chains else "(none)")

            res_remove = set()
            if het_sorted:
                sel = input("🧽 Select HET indices to remove for THIS file (1,2,...) or 'all' or blank to skip: ").strip().lower()
                if sel == "all":
                    res_remove = {name for name, _ in het_sorted}
                elif sel:
                    try:
                        chosen = {int(x.strip()) for x in sel.split(",") if x.strip()}
                        for i in chosen:
                            if 1 <= i <= len(het_sorted):
                                res_remove.add(het_sorted[i-1][0])
                            else:
                                print(f"  ⚠️ Ignoring out-of-range index: {i}")
                    except ValueError:
                        print("  ⚠️ Invalid index list, skipping automatic HET removal.")

            extra_res = input("➕ Extra 3-letter residue names to remove (optional): ").strip()
            if extra_res:
                res_remove |= {r.strip().upper() for r in extra_res.split(",") if r.strip()}

            if res_remove:
                print("🧼 Will remove HET residues (this file):", ", ".join(sorted(res_remove)))
            else:
                print("ℹ️ No HET residue names selected (this file).")

            chain_resp = input("🧩 Chains to remove for THIS file (comma-separated, blank to keep all): ").strip()
            chains_to_remove = {c.strip() for c in chain_resp.split(",") if c.strip()} if chain_resp else set()
            if chains_to_remove:
                print("✂️ Will remove chains (this file):", ", ".join(sorted(chains_to_remove)))
            else:
                print("ℹ️ Keeping all chains (this file).")

            per_file_choices[fname] = (res_remove, chains_to_remove)
            print()

    print(f"🔧 Using receptor preparer: {prep_kind}")

    # Process files
    success = 0
    fail = 0

    for fname in files:
        in_path = os.path.join(selected_folder, fname)
        base = os.path.splitext(fname)[0]
        cleaned_pdb = os.path.join(temp_cleaned_dir, base + ".clean.pdb")
        residues_to_remove, chains_to_remove = per_file_choices.get(fname, (set(), set()))

        try:
            clean_structure_to_pdb(in_path, cleaned_pdb, residues_to_remove, chains_to_remove)
            print(f"\n🧹 Cleaning {fname} → {os.path.basename(cleaned_pdb)}")
        except Exception as e:
            print(f"❌ Failed cleaning {fname}: {e}")
            fail += 1
            continue

        out_pdbqt = os.path.join(output_dir, base + ".converted.pdbqt")
        cmd = list(prep_base) + ["-r", cleaned_pdb, "-o", out_pdbqt]
        if prep_kind != "meeko":
            cmd += ["-A", "checkhydrogens"]
        env = prep_env if prep_env is not None else os.environ.copy()

        try:
            print(f"⚙️ Converting {fname} → {os.path.basename(out_pdbqt)}")
            subprocess.run(cmd, check=True, env=env)
            print(f"✅ Saved: {out_pdbqt}")
            success += 1
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed conversion on {fname}: {e}")
            fail += 1

    # Cleanup temp
    try:
        shutil.rmtree(temp_cleaned_dir)
    except Exception:
        pass

    # Final sweep: ensure ONLY .pdbqt remain in output folder
    for g in os.listdir(output_dir):
        if not g.lower().endswith(".pdbqt"):
            try:
                os.remove(os.path.join(output_dir, g))
            except Exception:
                pass

    print("\n🎉 Done.")
    print(f"   ✅ Converted: {success}")
    print(f"   ❌ Failed:    {fail}")
    print(f"📦 Output folder: {output_dir}")

if __name__ == "__main__":
    main()
