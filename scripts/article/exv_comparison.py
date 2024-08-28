import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as pe
import numpy as np
import sys
import os
from enum import Enum
from scipy.optimize import curve_fit

skip_similar_results = True

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
waxsis_folder = "output/waxsis/fitted"
match len(sys.argv):
    case 1: pass
    case 2:
        folder = sys.argv[1]
        if not os.path.exists(folder):
            print(f"Folder {sys.argv[2]} does not exist.")
            exit(1)
    case 3: 
        folder = sys.argv[1]
        waxsis_folder = "output/waxsis/"+sys.argv[2]
        if not os.path.exists(folder):
            print(f"Folder {sys.argv[1]} does not exist.")
            exit(1)
        if not os.path.exists(waxsis_folder):
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
    "HistogramManagerMTFFExplicit_voronoi:":        3, 
    "HistogramManagerMTFFExplicit_fitted_voronoi:": 4,
    "HistogramManagerMTFFExplicit_mf:":             5, 
    "HistogramManagerMTFFExplicit_fitted_mf:":      6,
    "HistogramManagerMTFFExplicit_vdw:":            7, 
    "HistogramManagerMTFFExplicit_fitted_vdw:":     8,
    "HistogramManagerMTFFExplicit_traube:":         9, 
    "HistogramManagerMTFFExplicit_fitted_traube:":  10,
    "HistogramManagerMTFFGrid:":                    11,
    "HistogramManagerMTFFGridSurface:":             12,
    "HistogramManagerMTFFGridSurface_fitted:":      13,
    "HistogramManagerMTFFGridSurface_fitted_300:":  14,
    "FoXS:":                                        15,
    "FoXS_fitted:":                                 16,
    "CRYSOL:":                                      17,
    "CRYSOL_fitted:":                               18,
    "Pepsi-SAXS:":                                  19,
    "Pepsi-SAXS_fitted:":                           20,
}

x_labels = [
    "Simple", 
    "Averaged", 
    "Averaged fitted", 
    "Explicit Voronoi", 
    "Voronoi", # Explicit Voronoi fitted
    "Explicit MF",
    "Gaussian spheres", # Explicit MF fitted
    "Explicit VDW",
    "van der Waals", # Explicit VDW fitted
    "Explicit Traube",
    "Traube", # Explicit Traube fitted
    "Grid", 
    "Grid surface",
    "Grid-based", 
    "Grid fitted 3.00", 
    "FoXS (ours)",
    "FoXS fitted (ours)",
    "CRYSOL (ours)",
    "CRYSOL fitted (ours)",
    "Pepsi-SAXS (ours)",
    "Pepsi-SAXS fitted (ours)",
    "CRYSOL", 
    "FoXS", 
    "Pepsi-SAXS",
    "WAXSiS"
]

class options(Enum):
    Simple = 0
    Average = 1
    Average_f = 2
    Explicit_VORONOI = 3
    Explicit_VORONOI_f = 4
    Explicit_MF = 5
    Explicit_MF_f = 6
    Explicit_VDW = 7
    Explicit_VDW_f = 8
    Explicit_TRAUBE = 9
    Explicit_TRAUBE_f = 10
    Grid = 11
    Grid_SURFACE = 12
    Grid_f_215 = 13
    Grid_f_300 = 14
    FoXS_o = 15
    FoXS_of = 16
    CRYSOL_o = 17
    CRYSOL_of = 18
    Pepsi_o = 19
    Pepsi_of = 20
    CRYSOL = 21
    FoXS = 22
    Pepsi = 23
    WAXSiS = 24

def waxsis_fit(folder):
    def chi2(ymodel):
        return np.sum(((data[:, 1] - ymodel) / data[:, 2]) ** 2)

    def perform_fit(fitdata, dof = 3):
        popt, _ = curve_fit(lambda x, a, b: a*x + b, fitdata[:, 1], data[:, 1], sigma=data[:, 2], absolute_sigma=True, p0=[1, 0])
        fitdata[:, 1] = popt[0] * fitdata[:, 1] + popt[1]
        chi2r = chi2(fitdata[:, 1]) / (len(data[:, 1]) - dof)
        return chi2r

    data_file = None
    fit_file = None
    for nested_folder in os.listdir(folder):
        fit_file = os.path.join(folder, nested_folder, "intensity.dat")
        for nested_file in os.listdir(os.path.join(folder, nested_folder, "envelope")):
            if nested_file[-4:] == ".dat":
                data_file = os.path.join(folder, nested_folder, "envelope", nested_file)
                continue

    if data_file is None or fit_file is None:
        print(f"waxsis_fit: Missing files in {folder}")
        return None

    data = np.loadtxt(data_file, skiprows=1)
    fitdata = np.loadtxt(fit_file, skiprows=0, comments="#", usecols=[0, 1])
    x = fitdata[:, 0]
    y = np.interp(data[:, 0], x, fitdata[:, 1])
    fitdata = np.vstack((data[:, 0], y)).T
    return perform_fit(fitdata)

size = []
foxs = []
pepsi = []
crysol = []
waxsis = []
y_labels = []   # file_name
data = []
ausaxs_length = 0
waxsis_folders = [os.path.join(waxsis_folder, folder) for folder in os.listdir(waxsis_folder) if os.path.isdir(os.path.join(waxsis_folder, folder))]
for i in range(len(options)):
    data.append([])
for file in os.listdir(folder):
    # open nested folders
    if os.path.isdir(os.path.join(folder, file)):
        found_waxsis = False
        found_ausaxs = False
        found_crysol = False
        found_pepsi = False
        found_foxs = False

        # search for matching waxsis data
        for waxsis_folder in waxsis_folders:
            if file in waxsis_folder:
                chi2r = waxsis_fit(waxsis_folder)
                if chi2r is not None:
                    waxsis.append(chi2r)
                else:
                    waxsis.append(0)
                found_waxsis = True
                break

        for nested_file in os.listdir(os.path.join(folder, file)):
            obj = os.path.join(folder, file, nested_file)
            if nested_file == "crysol.fit":
                with open(obj, "r") as f:
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
                found_crysol = True

            elif nested_file == "pepsi.fit":
                with open(obj, "r") as f:
                    tokens = f.readlines()[4].strip().split()
                    for i, token in enumerate(tokens):
                        if token.startswith("Chi"):
                            val = float(tokens[i+1])
                            if np.isnan(val) or val == None:
                                print(f"NaN in {file}")
                                exit(0)
                            pepsi.append(val)
                            break
                found_pepsi = True

            elif nested_file == "foxs.fit":
                with open(obj, "r") as f:
                    tokens = f.readlines()[1].strip().split()
                    for i, token in enumerate(tokens):
                        if token.startswith("Chi"):
                            val = float(tokens[i+2])
                            if np.isnan(val) or val == None:
                                print(f"NaN in {file}")
                                exit(0)
                            foxs.append(val)
                            break
                found_foxs = True

            elif nested_file == "ausaxs_chi2.txt":
                with open(obj, "r") as f:
                    tokens = f.readlines()
                    if ausaxs_length == 0:
                        ausaxs_length = len(tokens)
                    elif ausaxs_length != len(tokens):
                        print(f"Length mismatch in {file}")
                        exit(0)
                    for i, token in enumerate(tokens):
                        token = token.strip().split()
                        if token[0].startswith("size"):
                            size.append(int(token[1]))
                            continue

                        if token[0] in x_label_map:
                            val = float(token[1])
                            if np.isnan(val) or val == None:
                                print(f"NaN in {file}")
                                exit(0)
                            data[x_label_map[token[0]]].append(val)
                        else:
                            print(f"Unknown label {token[0]} in {file}")
                            exit(0)
                found_ausaxs = True

        name = file.split("_")[0]
        if name not in y_labels:
            y_labels.append(name)

        if not found_waxsis:
            waxsis.append(0)
        if not found_ausaxs:
            print(f"Missing AUSAXS data in {file}")
            exit(0)
        if not found_crysol:
            print(f"Missing CRYSOL data in {file}")
            exit(0)
        if not found_pepsi:
            print(f"Missing Pepsi data in {file}")
            exit(0)
        if not found_foxs:
            print(f"Missing Foxs data in {file}")
            exit(0)
        continue

def choose_data(choices: list[options]):
    data_out = []
    labels_out = []
    for choice in choices:
        if choice.value < options.CRYSOL.value:
            data_out.append(data[choice.value])
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
                case options.WAXSiS:
                    data_out.append(waxsis)
                    labels_out.append(x_labels[choice.value])
    return np.array(data_out).T, labels_out

###################################################
###                    SETUP                   ####
###################################################
size = np.array(size)
data_plot, x_labels_plot = choose_data([
    options.Simple, 
    # options.Average_f, 
    options.Explicit_MF_f,
    # options.Grid, 
    options.Grid_f_215, 
    # options.Grid_f_300, 
    options.CRYSOL, 
    options.FoXS, 
    options.Pepsi,
    options.WAXSiS
])

data_volumes, x_labels_volumes = choose_data([
    options.Explicit_VORONOI_f,
    options.Explicit_VDW_f,
    options.Explicit_TRAUBE_f
])
data_diff = data_volumes.copy()
x_labels_diff = x_labels_volumes
# sort by size
# indices = np.argsort(size)
# size = size[indices]
# data_pontius = data_pontius[indices]
# data_plot = data_plot[indices]
# y_labels = [y_labels[i] for i in indices]

# remove similar results
if skip_similar_results:
    indices = []
    reason_less_than_two = 0
    reason_similar = 0
    reason_too_large = 0
    max_chi2 = 20
    for i in range(len(data_plot)):
        chi2 = data_plot[i][0]
        add = False
        for j in range(1, len(data_plot[i])):
            if (max_chi2 < data_plot[i][j]) or data_plot[i][j] == 0:
                continue
            if (0.2 < abs(chi2 - data_plot[i][j])/chi2):
                add = True
                break
        for j in range(len(data_volumes[i])):
            if (max_chi2 < data_plot[i][j]) or data_plot[i][j] == 0:
                continue
            if (0.2 < abs(chi2 - data_plot[i][j])/chi2):
                add = True
                break
        if (data_plot[i] < 2).all() and (data_volumes[i] < 2).all():
            reason_less_than_two += 1
            add = False
        elif (max_chi2 < data_plot[i]).all() and (max_chi2 < data_volumes[i]).all():
            reason_too_large += 1
            add = False
        elif not add:
            reason_similar += 1

        if add:
            indices.append(i)

    print(f"Skipped {len(data_plot) - len(indices)} results")
    print(f"\t{reason_too_large} too large")
    print(f"\t{reason_less_than_two} less than 2")
    print(f"\t{reason_similar} too similar")
    data_plot = data_plot[indices]
    data_volumes = data_volumes[indices]
    y_labels = [y_labels[i] for i in indices]

tmp, _ = choose_data([options.Explicit_MF_f])
for i in range(len(data_diff)):
    data_diff[i] = tmp[i] - data_diff[i]
if skip_similar_results: data_diff = data_diff[indices]

bounds = [1, 1.5, 2, 3, 4, 5, 7.5, 10, 15, 25, 100]
cmap = plt.get_cmap('RdYlGn_r')
cmap.set_bad(color='white')
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

def round(x):
    if abs(x) < 5:
        return np.round(x, 2)
    if abs(x) < 10:
        return np.round(x, 1)
    return int(np.round(x))

def plot_one(data, x_labels, method):
    figsize_x = max(len(x_labels)+2, 8)
    figsize_y = max(len(data), 6)
    fig, ax = plt.subplots(figsize=(figsize_x, figsize_y))
    data[data[:, :] == 0] = np.NaN
    plt.imshow(data, interpolation='nearest', cmap=cmap, norm=norm)
    for i in range(len(y_labels)):
        for j in range(len(x_labels)):
            if np.isnan(data[i, j]):
                continue
            ax.text(j, i, round(data[i, j]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
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
            if np.isnan(data[j, i]):
                continue
            ax.text(j, i, round(data[j, i]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
        plt.axhline(i+0.5, color="black", linewidth=1)
    for i in range(len(y_labels)):
        plt.axvline(i+0.5, color="black", linewidth=1)

    plt.xticks(np.arange(len(y_labels)), y_labels, rotation=60, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(x_labels)), x_labels)
    plt.tight_layout()
    plt.savefig(f"output/fit_all_exv/{method}_comparison_flipped.png", dpi=300)

def plot_both():
    figsize_x = max(len(x_labels_plot)+len(x_labels_volumes)+2, 8)
    figsize_y = max(len(data_plot), 6)
    fig, ax = plt.subplots(1, 2, figsize=(figsize_x, figsize_y), sharey=True, gridspec_kw={'width_ratios': [len(x_labels_plot), len(x_labels_volumes) + 1]})

    plt.sca(ax[0])
    data_plot[data_plot[:, :] == 0] = np.NaN
    plt.imshow(data_plot, interpolation='nearest', cmap=cmap, norm=norm)
    for i in range(len(y_labels)):
        for j in range(len(x_labels_plot)):
            if np.isnan(data_plot[i, j]):
                continue
            ax[0].text(j, i, round(data_plot[i, j]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
        plt.axhline(i+0.5, color="black", linewidth=1)
    for i in range(len(x_labels_plot)):
        plt.axvline(i+0.5, color="black", linewidth=1)
    plt.xticks(np.arange(len(x_labels_plot)), x_labels_plot, rotation=60, ha="right", rotation_mode="anchor")
    ax[0].xaxis.set_ticks_position('none')
    ax[0].yaxis.set_ticks_position('none')

    plt.sca(ax[1])
    data_volumes[data_volumes[:, :] == 0] = np.NaN
    plt.imshow(data_volumes, interpolation='nearest', cmap=cmap, norm=norm)
    for i in range(len(y_labels)):
        for j in range(len(x_labels_volumes)):
            if np.isnan(data_volumes[i, j]):
                continue
            ax[1].text(j, i, round(data_volumes[i, j]), ha="center", va="center", color="w", fontsize=14, path_effects=[pe.withStroke(linewidth=1, foreground="black")])
        plt.axhline(i+0.5, color="black", linewidth=1)
    for i in range(len(x_labels_volumes)):
        plt.axvline(i+0.5, color="black", linewidth=1)
    plt.xticks(np.arange(len(x_labels_volumes)), x_labels_volumes, rotation=60, ha="right", rotation_mode="anchor")
    ax[1].xaxis.set_ticks_position('none')
    ax[1].yaxis.set_ticks_position('none')

    plt.sca(ax[0])
    plt.yticks(np.arange(len(y_labels)), y_labels)
    plt.tight_layout()
    plt.savefig(f"output/fit_all_exv/both.png", dpi=300)

plot_both()
if data_plot.shape[0] > 0:
    plot_one(data_plot, x_labels_plot, "Main")

if data_volumes.shape[0] > 0:
    plot_one(data_volumes, x_labels_volumes, "Other volumes")

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

plt.title("Volume differences")
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

plt.title("Volume differences")
plt.xticks(np.arange(len(y_labels)), y_labels, rotation=60, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(x_labels_diff)), x_labels_diff)
plt.tight_layout()
plt.savefig(f"output/fit_all_exv/difference_flipped.png", dpi=300)
