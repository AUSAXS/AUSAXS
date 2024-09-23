from plot_extensions import load_file 
from plot_helper import *

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter

params = {
    'legend.fontsize': 24,
    'figure.figsize': (10, 8),
    'axes.labelsize': 24,
    'axes.titlesize': 24,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'lines.markersize': 10,
    'lines.linewidth': 3
}
plt.rcParams.update(params)

# fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
# plt.sca(ax[0])
# data = load_file("output/vary_grid_cell_size/intensities.png.plot")
# for d in data[:4]:
#     if isinstance(d, Dataset):
#         plot_dataset(d)

# ax[0].set_xscale("linear")
# ax[0].set_ylim(1e-4, 1e-2)

# plt.sca(ax[1])
# data = load_file("output/vary_grid_cell_size/difference.png.plot")
# for d in data[:4]:
#     if isinstance(d, Dataset):
#         plot_dataset(d)

# ax[1].set_xscale("linear")
# ax[1].set_ylim(0.5, 1.05)

plt.figure(figsize=(10, 6))
data = load_file("output/vary_grid_cell_size/difference.png.plot")
for d in data[:4]:
    if isinstance(d, Dataset):
        plot_dataset(d)

plt.legend(["$d_{cell}$ = 1.0 Å", "$d_{cell}$ = 1.5 Å", "$d_{cell}$ = 2.0 Å", "$d_{cell}$ = 2.5 Å"])
plt.axhline(1, color="black", linestyle="-")
plt.xscale("linear")
plt.yscale("log")
plt.ylabel("Baseline ratio")
plt.xlabel("$q\ [\AA^{-1}]$")
plt.tight_layout()
plt.gca().yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# plt.ylim(0.5, 1.05)
plt.savefig("output/vary_grid_cell_size/difference.png", dpi=300)
plt.show()