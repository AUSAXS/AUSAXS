from plot_extensions import load_file 
from plot_helper import *

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter

title_font_size = 20
params = {
    'legend.fontsize': 16,
    'figure.figsize': (10, 8),
    'axes.labelsize': 20,
    'axes.titlesize': 16,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'lines.markersize': 10,
    'lines.linewidth': 3
}
plt.rcParams.update(params)

fig, ax = plt.subplots(2, 2, figsize=(16, 10), sharex=True, sharey=True, gridspec_kw={'hspace': 0.075, 'wspace': 0.05})
data = load_file("output/bin_size_analysis/unweighted_unstructured.png.plot")
plt.sca(ax[0, 0])
for d in data:
    if isinstance(d, Dataset):
        plot_dataset(d)
plt.xscale("linear")
plt.yscale("log")
plt.xlabel("")
plt.ylabel("Relative deviation")
plt.yticks([0.5, 1, 2], ["0.5", "1", "2"])
plt.minorticks_off()
plt.legend(["$w=0.50$", "$w=0.25$", "$w=0.15$", "$w=0.10$", "$w=0.05$"], loc='upper right', ncol=3)

plt.sca(ax[0, 1])
data = load_file("output/bin_size_analysis/unweighted_structured.png.plot")
for d in data:
    if isinstance(d, Dataset):
        plot_dataset(d)
plt.xscale("linear")
plt.yscale("log")
plt.xlabel("")
plt.ylabel("")
plt.yticks([0.5, 1, 2], ["0.5", "1", "2"])
plt.minorticks_off()
plt.gca().get_legend().remove()

plt.sca(ax[1, 0])
data = load_file("output/bin_size_analysis/weighted_unstructured.png.plot")
for d in data:
    if isinstance(d, Dataset):
        plot_dataset(d)
plt.xscale("linear")
plt.yscale("log")
plt.xlabel("q [$\AA^{-1}$]")
plt.ylabel("Relative deviation")
plt.yticks([0.5, 1, 2], ["0.5", "1", "2"])
plt.minorticks_off()
plt.gca().get_legend().remove()

plt.sca(ax[1, 1])
data = load_file("output/bin_size_analysis/weighted_structured.png.plot")
for d in data:
    if isinstance(d, Dataset):
        plot_dataset(d)
plt.xscale("linear")
plt.yscale("log")
plt.xlabel("q [$\AA^{-1}$]")
plt.ylabel("")
plt.yticks([0.5, 1, 2], ["0.5", "1", "2"])
plt.minorticks_off()
plt.gca().get_legend().remove()

plt.text(1.05, 0.5, "Unweighted", fontsize=title_font_size, ha='center', va='center', transform=ax[0, 1].transAxes, rotation=270)
plt.text(1.05, 0.5, "Weighted", fontsize=title_font_size, ha='center', va='center', transform=ax[1, 1].transAxes, rotation=270)
plt.text(0.5, 1.05, "Structured", fontsize=title_font_size, ha='center', va='center', transform=ax[0, 1].transAxes)
plt.text(0.5, 1.05, "Unstructured", fontsize=title_font_size, ha='center', va='center', transform=ax[0, 0].transAxes)

# plt.tight_layout()
plt.savefig("output/bin_size_analysis/weighted_bins.png", dpi=300)
plt.show()