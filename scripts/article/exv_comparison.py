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
match len(sys.argv):
    case 1:
        print("missing one required argument: folder")
        exit(1)

    case 2:
        folder = sys.argv[1]
        if not os.path.exists(folder):
            print(f"Folder {sys.argv[1]} does not exist.")
            exit(1)
        
if folder is None:
    print("Folder not found.")
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
x_labels = ["Simple", "Averaged", "Averaged & fitted", "Explicit", "Explicit & fitted", "Grid"]   # method
y_labels = []   # file_name
for file in os.listdir(folder):
    if not file.endswith(".txt"):
        continue

    # filenames are TRAUBE_name.txt or PONTIUS_name.txt; we want to keep only the name
    y_labels.append(file.split("_")[1].split(".")[0])
    with open(os.path.join(folder, file), "r") as f:
        chi2s = np.array([0.0]*6)
        for line in f:
            line = line.strip().split()
            chi2s[x_label_map[line[0]]] = float(line[1])
        data.append(chi2s)

data = np.array(data)
fig, ax = plt.subplots()

bounds = [1, 1.5, 2, 3, 4, 5, 7.5, 10, 15, 25, 100]
cmap = plt.get_cmap('RdYlGn_r')
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

def round(x):
    if x < 5:
        return np.round(x, 2)
    if x < 10:
        return np.round(x, 1)
    return int(np.round(x))

plt.imshow(data, interpolation='nearest', cmap=cmap, norm=norm)
for i in range(len(y_labels)):
    for j in range(len(x_labels)):
        text = ax.text(j, i, round(data[i, j]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
    plt.axhline(i+0.5, color="black", linewidth=1)
for i in range(len(x_labels)):
    plt.axvline(i+0.5, color="black", linewidth=1)

plt.xticks(np.arange(len(x_labels)), x_labels, rotation=60, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(y_labels)), y_labels)
plt.tight_layout()
plt.show()