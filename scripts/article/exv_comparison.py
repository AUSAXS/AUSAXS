import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as pe
import numpy as np
import sys
import os

params = {
    'legend.fontsize': 14,
    'figure.figsize': (10, 8),
    'axes.labelsize': 14,
    'axes.titlesize':14,
    'xtick.labelsize':11,
    'ytick.labelsize':11,
    'lines.markersize': 5
}

folder = ""
method = ""
match len(sys.argv):
    case 3:
        folder = sys.argv[1]
        if not os.path.exists(folder):
            print(f"Folder {sys.argv[2]} does not exist.")
            exit(1)
        method = sys.argv[2]
        if method != "TRAUBE" and method != "PONTIUS":
            print("Method must be either TRAUBE or PONTIUS.")
            exit(1)

    case _:
        print("Usage: python exv_comparison.py <folder> <method>")
        exit(1)

# iterate through all files in the directory
x_label_map = {
    "HistogramManagerMT:": 0, 
    "HistogramManagerMTFFAvg:": 1, 
    "HistogramManagerMTFFAvg_fitted:": 2, 
    "HistogramManagerMTFFExplicit:": 3, 
    "HistogramManagerMTFFExplicit_fitted:": 4,
    "HistogramManagerMTFFGrid:": 5
}
data = []       # indexing: [file_name, method]
foxs = []
pepsi = []
crysol = []
x_labels = ["Simple", "Averaged", "Averaged & fitted", "Explicit", "Explicit & fitted", "Grid"]   # method
y_labels = []   # file_name
for file in os.listdir(folder):
    # open nested folders
    if os.path.isdir(os.path.join(folder, file)):
        for nested_file in os.listdir(os.path.join(folder, file)):
            if not (nested_file.endswith(".txt") or nested_file.endswith("fit")):
                continue
            
            if nested_file == "crysol.fit":
                with open(os.path.join(folder, file, nested_file), "r") as f:
                    tokens = f.readline().strip().split()
                    for i, token in enumerate(tokens):
                        if token.startswith("Chi"):
                            if (len(tokens)) == i+2:
                                if tokens[i+1].startswith("*"):
                                    print(f"Error in {file}: {tokens[i+1]}")
                                crysol.append(float(tokens[i+1]))
                            else:
                                if token.split(":")[1].startswith("*"):
                                    print(f"Error in {file}: {token.split(':')[1]}")
                                crysol.append(float(token.split(":")[1]))
                            break

            elif nested_file == "pepsi.fit":
                with open(os.path.join(folder, file, nested_file), "r") as f:
                    tokens = f.readlines()[4].strip().split()
                    for i, token in enumerate(tokens):
                        if token.startswith("Chi"):
                            pepsi.append(float(tokens[i+1]))
                            break

            elif nested_file == "foxs.fit":
                with open(os.path.join(folder, file, nested_file), "r") as f:
                    tokens = f.readlines()[1].strip().split()
                    for i, token in enumerate(tokens):
                        if token.startswith("Chi"):
                            foxs.append(float(tokens[i+2]))
                            break
        
            if not nested_file.startswith(method):
                continue

            # filenames are TRAUBE_name.txt or PONTIUS_name.txt; we want to keep only the name
            y_labels.append(nested_file.split("_")[1].split(".")[0])
            chi2s = np.array([0.0]*6)
            with open(os.path.join(folder, file, nested_file), "r") as f:
                for line in f:
                    line = line.strip().split()
                    chi2s[x_label_map[line[0]]] = float(line[1])
                data.append(chi2s)
        continue

data = np.array(data)
if method == "TRAUBE":
    data = np.concatenate((data, np.array([crysol, foxs, pepsi]).T), axis=1)
    x_labels.extend(["CRYSOL", "FoXS", "Pepsi-SAXS"])
bounds = [1, 1.5, 2, 3, 4, 5, 7.5, 10, 15, 25, 100]
cmap = plt.get_cmap('RdYlGn_r')
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

def round(x):
    if x < 5:
        return np.round(x, 2)
    if x < 10:
        return np.round(x, 1)
    return int(np.round(x))

fig, ax = plt.subplots(figsize=(len(x_labels), len(data)))
plt.imshow(data, interpolation='nearest', cmap=cmap, norm=norm)
for i in range(len(y_labels)):
    for j in range(len(x_labels)):
        text = ax.text(j, i, round(data[i, j]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
    plt.axhline(i+0.5, color="black", linewidth=1)
for i in range(len(x_labels)):
    plt.axvline(i+0.5, color="black", linewidth=1)

plt.title(f"Comparison with {method} volume")
plt.xticks(np.arange(len(x_labels)), x_labels, rotation=60, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(y_labels)), y_labels)
plt.tight_layout()
plt.savefig(f"output/fit_all_exv/{method}_comparison.png", dpi=300)

# flipped axes
fig, ax = plt.subplots(figsize=(len(data), len(x_labels)))
plt.imshow(data.T, interpolation='nearest', cmap=cmap, norm=norm)
for i in range(len(x_labels)):
    for j in range(len(y_labels)):
        text = ax.text(j, i, round(data[j, i]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
    plt.axhline(i+0.5, color="black", linewidth=1)
for i in range(len(y_labels)):
    plt.axvline(i+0.5, color="black", linewidth=1)

plt.title(f"Comparison with {method} volume")
plt.xticks(np.arange(len(y_labels)), y_labels, rotation=60, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(x_labels)), x_labels)
plt.tight_layout()
plt.savefig(f"output/fit_all_exv/{method}_comparison_flipped.png", dpi=300)
plt.show()