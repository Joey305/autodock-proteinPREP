# import os
# import csv

# results_dir = os.path.join(os.getcwd(), "Docking_Results")
# output_csv = "vina_docking_scores.csv"

# rows = [("Receptor", "Ligand", "Pose", "Binding_Affinity_kcal_per_mol")]

# print(f"📂 Scanning: {results_dir}")

# for receptor in os.listdir(results_dir):
#     receptor_path = os.path.join(results_dir, receptor)
#     if not os.path.isdir(receptor_path):
#         continue

#     for ligand in os.listdir(receptor_path):
#         ligand_path = os.path.join(receptor_path, ligand)
#         log_file = os.path.join(ligand_path, "log.txt")

#         if not os.path.isfile(log_file):
#             print(f"❌ Missing log file: {log_file}")
#             continue

#         print(f"📄 Parsing {receptor} ← {ligand}")

#         with open(log_file, "r") as f:
#             lines = f.readlines()

#         in_results = False
#         for line in lines:
#             if line.strip().startswith("mode |"):
#                 in_results = True
#                 continue
#             if in_results:
#                 if not line.strip() or line.startswith("-"):
#                     continue
#                 parts = line.strip().split()
#                 if len(parts) >= 2 and parts[0].isdigit():
#                     pose = parts[0]
#                     score = parts[1]
#                     rows.append((receptor, ligand, pose, score))

# # Save results
# with open(output_csv, "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerows(rows)

# print(f"✅ CSV created: {output_csv}")

import os
import csv

results_dir = os.path.join(os.getcwd(), "Docking_Results")
output_csv = "vina_docking_scores_expanded.csv"

rows = [("Receptor", "Ligand", "SMILE_POSE", "Pose", "Binding_Affinity_kcal_per_mol")]

print(f"📂 Scanning: {results_dir}")

for receptor in os.listdir(results_dir):
    receptor_path = os.path.join(results_dir, receptor)
    if not os.path.isdir(receptor_path):
        continue

    for ligand in os.listdir(receptor_path):
        ligand_path = os.path.join(receptor_path, ligand)
        log_file = os.path.join(ligand_path, "log.txt")

        if not os.path.isfile(log_file):
            print(f"❌ Missing log file: {log_file}")
            continue

        print(f"📄 Parsing {receptor} ← {ligand}")

        try:
            lig_prefix, smile_pose = ligand.split("__", 1)
        except ValueError:
            print(f"⚠️ Could not split ligand name: {ligand}")
            lig_prefix = ligand
            smile_pose = "unknown"

        with open(log_file, "r") as f:
            lines = f.readlines()

        in_results = False
        for line in lines:
            if line.strip().startswith("mode |"):
                in_results = True
                continue
            if in_results:
                if not line.strip() or line.startswith("-"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0].isdigit():
                    pose = parts[0]
                    score = parts[1]
                    rows.append((receptor, lig_prefix, smile_pose, pose, score))

# Save results
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(rows)

print(f"✅ CSV created: {output_csv}")
