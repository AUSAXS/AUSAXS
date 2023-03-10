import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.optimize import curve_fit
import argparse

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
files = [os.path.join(folder, f) for f in files]

# plot each file in log-log scale, staggered by 10
for i, f in enumerate(files):
    # get the stem of f without the path or extension
    stem = os.path.splitext(os.path.basename(f))[0]
    print("\tParsing file: " + stem)

    # load the data
    data = np.loadtxt(f, skiprows=0, comments=["@", "#", "&"], usecols=[0, 1])
    data[:, 1] = data[:, 1] * 10**i

    # plot the data
    plt.plot(data[:, 0], data[:, 1], label=stem.lower(), alpha=0.5)

# add a legend and save the plot
plt.legend()
plt.yscale("log")
plt.ylabel("Intensity [staggered]")
plt.savefig(os.path.join(folder, "log.png"))
plt.xscale("log")
plt.xlabel("q [Å$^{-1}$]")
plt.savefig(os.path.join(folder, "loglog.png"))
plt.show()