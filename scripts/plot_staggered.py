import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.optimize import curve_fit
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

# plot each file in log-log scale, staggered by 10
skip = 5
mainfig = plt.figure() # main figure with all plots
subfig = plt.figure()  # subfigure with only skipfig at a time
skipfig = plt.figure() # figure plotting only every skipfig file
counter = 0
for i, f in enumerate(files, start=1):
    # get the stem of f without the path or extension
    stem = os.path.splitext(os.path.basename(f))[0]
    print("\tParsing file: " + stem)

    # load the data
    data = np.loadtxt(f, skiprows=0, comments=["@", "#", "&"], usecols=[0, 1])
    # data[:, 1] = data[:, 1] * 2**i

    # remove 10 last points
    data = np.delete(data, np.s_[-10:], axis=0)

    # plot the data
    plt.figure(mainfig.number)
    plt.plot(data[:, 0], data[:, 1], label=stem.lower(), alpha=0.5)
    plt.figure(subfig.number)
    plt.plot(data[:, 0], data[:, 1], label=stem.lower(), alpha=0.5)
    if i % skip == 0:
        counter += 1
        # save subfig
        plt.legend()
        plt.yscale("log")
        plt.ylabel("Intensity [staggered]")
        plt.savefig(os.path.join(folder, "log.png"))
        plt.xscale("log")
        plt.xlabel("q [Å$^{-1}$]")
        plt.savefig(os.path.join(folder, f"loglog_{counter}.png"))
        plt.close(subfig)
        subfig = plt.figure()
        
        plt.figure(skipfig.number)
        plt.plot(data[:, 0], data[:, 1], label=stem.lower(), alpha=0.5)

# add a legend and save the plot
plt.figure(mainfig.number)
plt.legend()
plt.yscale("log")
plt.ylabel("Intensity [staggered]")
plt.savefig(os.path.join(folder, "log.png"))
plt.xscale("log")
plt.xlabel("q [Å$^{-1}$]")
plt.savefig(os.path.join(folder, "loglog.png"))

plt.figure(skipfig.number)
plt.legend()
plt.yscale("log")
plt.ylabel("Intensity [staggered]")
plt.savefig(os.path.join(folder, "log_skip.png"))
plt.xscale("log")
plt.xlabel("q [Å$^{-1}$]")
plt.savefig(os.path.join(folder, "loglog_skip.png"))