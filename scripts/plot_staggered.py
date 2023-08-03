import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import argparse

params = {
    'legend.fontsize': 14,
    'figure.figsize': (20, 14),
    'axes.labelsize': 14,
    'axes.titlesize':14,
    'xtick.labelsize':11,
    'ytick.labelsize':11
}
plt.rcParams.update(params)

parser = argparse.ArgumentParser(description="Plot data files.")
parser.add_argument("folder",          help="The folder containing the data files.")
parser.add_argument("-o", "--output",  help="The output path.")
args = parser.parse_args()

folder: str = args.folder
if not os.path.exists(folder):
    print(f"Folder {folder} does not exist.")
    sys.exit(1)

if args.output:
    folder = args.output
    if not os.path.exists(folder):
        os.mkdir(folder)

# get all the files in the folder
files = [f for f in os.listdir(folder) if f.endswith(('.fit', ".xvg"))]
files = sorted([os.path.join(folder, f) for f in files], key=lambda x: int(x.split("/")[-1].split(".")[0]))

# define data struct
class Data:
    def __init__(self, name, x, y, dy):
        self.x = x
        self.y = y
        self.dy = dy
        self.name = name

    def len(self):
        return len(self.x)

# load data
data : list[Data] = []
for f in files:
    # get the stem of f without the path or extension
    print("\tParsing file: " + f)
    stem = os.path.splitext(os.path.basename(f))[0]

    # load the data
    tmp = np.loadtxt(f, skiprows=0, comments=["@", "#", "&"], usecols=[0, 1, 2])
    data.append(Data(stem, tmp[:, 0], tmp[:, 1], tmp[:, 2]))

# calculate the weighted mean of the y values
avg = Data("average", np.zeros(len(data[0].x)), np.zeros(len(data[0].x)), np.zeros(len(data[0].x)))
for i in range(data[0].len()):
    w = 0
    yw = 0
    for d in data:
        w += 1 / d.dy[i]**2
        yw += d.y[i] / d.dy[i]**2

    avg.x[i] = data[0].x[i]
    avg.y[i] = yw / w
    avg.dy[i] = np.sqrt(1 / w) * np.var([d.y[i] for d in data]) # standard error of the weighted mean

########################################################
###                        ALL                       ###
########################################################
plt.figure()
for d in data:
    plt.plot(d.x, d.y, c="k", alpha=0.4)
plt.plot(avg.x, avg.y, c="r", label="average")
plt.xlabel("q [Å$^{-1}$]")
plt.ylabel("Intensity [arb]")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig(os.path.join(folder, "all.png"))

########################################################
###                     SINGLE                       ###
########################################################
chi2s = []
for i, _ in enumerate(data):
    # calculate the autocorrelation with the other datasets
    autocorr = []
    for j in range(i, len(data)):
        auto = np.correlate(data[i].y, data[j].y, mode="full")

    # calculate chi2 with the average
    chi2 = np.sum((data[i].y - avg.y)**2 / data[i].dy**2)
    chi2s.append(chi2)

    # plot the fit in an upper subplot, and the residuals in a lower subplot
    fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    ax[0].errorbar(data[i].x, data[i].y, yerr=data[i].dy, c="k", label="data")
    ax[0].plot(avg.x, avg.y, c="r", label="average")
    ax[0].set_ylabel("Intensity [arb]")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].legend()

    # plot the residuals
    ax[1].axhline(0, c="k", lw=0.5)
    ax[1].plot(data[i].x, (data[i].y - avg.y) / data[i].dy, ".k", label=f"$\chi^2$ = {chi2:.3f}")
    ax[1].set_xlabel("q [Å$^{-1}$]")
    ax[1].set_ylabel("Residuals")
    ax[1].set_xscale("log")
    ax[1].legend()

    # save the figure
    plt.savefig(os.path.join(folder, f"conv_{data[i].name}.png"))
    plt.close()

# convergence plot
plt.figure()
x = np.arange(0, len(chi2s) * 0.5, 0.5)
plt.plot(x, chi2s, ".k")
plt.xlabel("Simulation time [ns]")
plt.ylabel("$\chi^2$")
plt.savefig(os.path.join(folder, "convergence.png"))

########################################################
###                 AUTOCORRELATION                  ###
########################################################
# calculate the autocorrelation of the y values
autocorr = Data("autocorrelation", np.zeros(len(data[0].x)), np.zeros(len(data[0].x)), np.zeros(len(data[0].x)))
for i in range(data[0].len()):
    w = 0
    yw = 0
    for d in data:
        w += 1 / d.dy[i]**2
        yw += d.y[i] / d.dy[i]**2

    autocorr.x[i] = data[0].x[i]
    autocorr.y[i] = yw / w
    autocorr.dy[i] = np.sqrt(1 / w) * np.var([d.y[i] for d in data]) # standard error of the weighted mean

# plot the autocorrelation
plt.figure()
plt.plot(autocorr.x, autocorr.y, c="k")
plt.xlabel("q [Å$^{-1}$]")
plt.ylabel("Autocorrelation")
plt.xscale("log")
plt.savefig(os.path.join(folder, "autocorrelation.png"))

# plot each file in log-log scale, staggered by 10
# skip = 5
# mainfig = plt.figure() # main figure with all plots
# subfig = plt.figure()  # subfigure with only skipfig at a time
# skipfig = plt.figure() # figure plotting only every skipfig file
# counter = 0
# for i, f in enumerate(files, start=1):
#     # get the stem of f without the path or extension
#     stem = os.path.splitext(os.path.basename(f))[0]
#     print("\tParsing file: " + stem)

#     # load the data
#     data = np.loadtxt(f, skiprows=0, comments=["@", "#", "&"], usecols=[0, 1])
#     # data[:, 1] = data[:, 1] * 2**i

#     # remove 10 last points
#     # data = np.delete(data, np.s_[10:-10], axis=0)

#     # plot the data
#     plt.figure(mainfig.number)
#     plt.plot(data[:, 0], data[:, 1], c="k", alpha=0.1)
#     if (i == 2):
#         plt.plot(data[:, 0], data[:, 1], label=stem.lower(), c="r", alpha=0.5)

#     plt.figure(subfig.number)
#     plt.plot(data[:, 0], data[:, 1], label=stem.lower(), alpha=0.5)
#     if i % skip == 0:
#         counter += 1
#         # save subfig
#         plt.legend()
#         plt.yscale("log")
#         plt.ylabel("Intensity [staggered]")
#         plt.savefig(os.path.join(folder, "log.png"))
#         plt.xscale("log")
#         plt.xlabel("q [Å$^{-1}$]")
#         plt.savefig(os.path.join(folder, f"loglog_{counter}.png"))
#         plt.close(subfig)
#         subfig = plt.figure()

#         plt.figure(skipfig.number)
#         plt.plot(data[:, 0], data[:, 1], label=stem.lower(), alpha=0.5)

# # add a legend and save the plot
# plt.figure(mainfig.number)
# plt.legend()
# plt.yscale("log")
# plt.ylabel("Intensity [staggered]")
# plt.savefig(os.path.join(folder, "log.png"))
# plt.xscale("log")
# plt.xlabel("q [Å$^{-1}$]")
# plt.savefig(os.path.join(folder, "loglog.png"))

# plt.figure(skipfig.number)
# plt.legend()
# plt.yscale("log")
# plt.ylabel("Intensity [staggered]")
# plt.savefig(os.path.join(folder, "log_skip.png"))
# plt.xscale("log")
# plt.xlabel("q [Å$^{-1}$]")
# plt.savefig(os.path.join(folder, "loglog_skip.png"))