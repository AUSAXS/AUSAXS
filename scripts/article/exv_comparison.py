import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as pe
import numpy as np
import sys
import os
from enum import Enum

params = {
    'legend.fontsize': 14,
    'figure.figsize': (10, 8),
    'axes.labelsize': 14,
    'axes.titlesize':14,
    'xtick.labelsize':11,
    'ytick.labelsize':11,
    'lines.markersize': 5
}

folder = "output/fit_all_exv/"
match len(sys.argv):
    case 1: pass
    case 2:
        folder = sys.argv[1]
        if not os.path.exists(folder):
            print(f"Folder {sys.argv[2]} does not exist.")
            exit(1)

    case _:
        print("Usage: python exv_comparison.py <folder>")
        exit(1)

# iterate through all files in the directory
x_label_map = {
    "HistogramManagerMT:":                          0, 
    "HistogramManagerMTFFAvg:":                     1, 
    "HistogramManagerMTFFAvg_fitted:":              2, 
    "HistogramManagerMTFFExplicit:":                3, 
    "HistogramManagerMTFFExplicit_fitted:":         4,
    "HistogramManagerMTFFGrid:":                    5,
    "HistogramManagerMTFFGridSurface:":             6,
    "HistogramManagerMTFFGridSurface_fitted_215:":  7,
    "HistogramManagerMTFFGridSurface_fitted_300:":  8,
    "FoXS:":                                        9,
    "FoXS_fitted:":                                10,
    "CRYSOL:":                                     11,
    "CRYSOL_fitted:":                              12,
    "Pepsi-SAXS:":                                 13,
    "Pepsi-SAXS_fitted:":                          14,
}

x_labels = [
    "Simple", 
    "Averaged", 
    "Averaged fitted", 
    "Explicit", 
    "Explicit fitted", 
    "Grid", 
    "Grid surface",
    "Grid fitted 2.15", 
    "Grid fitted 3.00", 
    "FoXS (ours)",
    "FoXS fitted (ours)",
    "CRYSOL (ours)",
    "CRYSOL fitted (ours)",
    "Pepsi-SAXS (ours)",
    "Pepsi-SAXS fitted (ours)",
    "CRYSOL", 
    "FoXS", 
    "Pepsi-SAXS"
]

class options(Enum):
    Simple = 0
    Average = 1
    Average_f = 2
    Explicit = 3
    Explicit_f = 4
    Grid = 5
    Grid_SURFACE = 6
    Grid_f_215 = 7
    Grid_f_300 = 8
    FoXS_o = 9
    FoXS_of = 10
    CRYSOL_o = 11
    CRYSOL_of = 12
    Pepsi_o = 13
    Pepsi_of = 14
    CRYSOL = 15
    FoXS = 16
    Pepsi = 17

data = {"TRAUBE": [], "PONTIUS": []}       # indexing: [file_name, method]
size = []
foxs = []
pepsi = []
crysol = []
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
        
            if not (nested_file.startswith("TRAUBE") or nested_file.startswith("PONTIUS")):
                continue

            # filenames are TRAUBE_name.txt or PONTIUS_name.txt; we want to keep only the name
            tokens = nested_file.split("_")
            name = tokens[1].split(".")[0]
            if (name not in y_labels):
                y_labels.append(name)
            chi2s = np.array([0.0]*len(x_label_map))
            with open(os.path.join(folder, file, nested_file), "r") as f:
                for line in f:
                    line = line.strip().split()
                    if (line[0].startswith("size")):
                        if (int(line[1]) not in size):
                            size.append(int(line[1]))
                    else:
                        if (line[0] in x_label_map and x_label_map[line[0]] is not None):
                            chi2s[x_label_map[line[0]]] = float(line[1])
                        else:
                            print(f"Unknown method: {line[0]}")
                data[tokens[0]].append(chi2s)
        continue

def choose_data(method: str, choices: list[options]):
    if method != "TRAUBE" and method != "PONTIUS":
        print(f"Unknown method: {method}")
        exit(1)

    data_out = []
    labels_out = []
    method_data = np.array(data[method])
    if len(method_data) == 0: return np.array(data_out), labels_out
    for choice in choices:
        if choice.value < options.CRYSOL.value:
            data_out.append(method_data[:, choice.value])
            labels_out.append(x_labels[choice.value])
        else:
            match choice:
                case options.CRYSOL:
                    data_out.append(crysol)
                    labels_out.append(x_labels[choice.value])
                case options.FoXS:
                    data_out.append(foxs)
                    labels_out.append(x_labels[choice.value])
                case options.Pepsi:
                    data_out.append(pepsi)
                    labels_out.append(x_labels[choice.value])
    return np.array(data_out).T, labels_out

###################################################
###                    SETUP                   ####
###################################################
size = np.array(size)
data_pontius = np.array(data["PONTIUS"])
data_traube = np.array(data["TRAUBE"])
data_diff = None

data_traube, x_labels_traube = choose_data("TRAUBE", [
    options.Simple, 
    options.Average_f, 
    options.Explicit_f, 
    options.Grid, 
    options.Grid_f_215, 
    options.Grid_f_300, 
    options.CRYSOL, 
    options.FoXS, 
    options.Pepsi
])

data_pontius, x_labels_pontius = choose_data("PONTIUS", [
    options.Average_f, 
    options.Explicit_f
])

# sort by size
# indices = np.argsort(size)
# size = size[indices]
# data_pontius = data_pontius[indices]
# data_traube = data_traube[indices]
# y_labels = [y_labels[i] for i in indices]

data_diff, x_labels_diff = choose_data("TRAUBE", [
    options.Average_f, 
    options.Explicit_f
])
data_diff -= data_pontius

bounds = [1, 1.5, 2, 3, 4, 5, 7.5, 10, 15, 25, 100]
cmap = plt.get_cmap('RdYlGn_r')
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

def round(x):
    if x < 5:
        return np.round(x, 2)
    if x < 10:
        return np.round(x, 1)
    return int(np.round(x))

def plot_one(data, x_labels, method):
    figsize_x = max(len(x_labels)+2, 8)
    figsize_y = max(len(data), 6)
    fig, ax = plt.subplots(figsize=(figsize_x, figsize_y))
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
    fig, ax = plt.subplots(figsize=(figsize_y, figsize_x))
    plt.imshow(data.T, interpolation='nearest', cmap=cmap, norm=norm)
    for i in range(len(x_labels)):
        for j in range(len(y_labels)):
            ax.text(j, i, round(data[j, i]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
        plt.axhline(i+0.5, color="black", linewidth=1)
    for i in range(len(y_labels)):
        plt.axvline(i+0.5, color="black", linewidth=1)

    plt.title(f"Comparison with {method} volume")
    plt.xticks(np.arange(len(y_labels)), y_labels, rotation=60, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(x_labels)), x_labels)
    plt.tight_layout()
    plt.savefig(f"output/fit_all_exv/{method}_comparison_flipped.png", dpi=300)

if data_traube.shape[0] > 0:
    plot_one(data_traube, x_labels_traube, "Traube")

if data_pontius.shape[0] > 0:
    plot_one(data_pontius, x_labels_pontius, "Pontius")


###################################################
###                 DIFFERENCE                 ####
###################################################
if data_diff is None:
    exit(0)

figsize_x = max(len(x_labels_diff)+2, 8)
figsize_y = max(len(data_diff), 6)

bounds = [-10, -5, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 5, 10]
cmap = plt.get_cmap('RdYlGn_r')
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(figsize_x, figsize_y))
plt.imshow(data_diff, interpolation='nearest', cmap=cmap, norm=norm)
for i in range(len(y_labels)):
    for j in range(len(x_labels_diff)):
        ax.text(j, i, round(data_diff[i, j]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
    plt.axhline(i+0.5, color="black", linewidth=1)
for i in range(len(x_labels_diff)):
    plt.axvline(i+0.5, color="black", linewidth=1)

plt.title("Difference between Traube and Pontius")
plt.xticks(np.arange(len(x_labels_diff)), x_labels_diff, rotation=60, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(y_labels)), y_labels)
plt.tight_layout()
plt.savefig(f"output/fit_all_exv/difference.png", dpi=300)

# flipped axes
fig, ax = plt.subplots(figsize=(figsize_y, figsize_x))
plt.imshow(data_diff.T, interpolation='nearest', cmap=cmap, norm=norm)
for i in range(len(x_labels_diff)):
    for j in range(len(y_labels)):
        ax.text(j, i, round(data_diff[j, i]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
    plt.axhline(i+0.5, color="black", linewidth=1)
for i in range(len(y_labels)):
    plt.axvline(i+0.5, color="black", linewidth=1)

plt.title("Difference between Traube and Pontius")
plt.xticks(np.arange(len(y_labels)), y_labels, rotation=60, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(x_labels_diff)), x_labels_diff)
plt.tight_layout()
plt.savefig(f"output/fit_all_exv/difference_flipped.png", dpi=300)
