import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.optimize import curve_fit
import argparse

elinewidth = 2
params = {
    'legend.fontsize': 20,
    'figure.figsize': (10, 8),
    'axes.labelsize': 20,
    'axes.titlesize':20,
    'xtick.labelsize':14,
    'ytick.labelsize':14,
    'lines.linewidth': 2,
    'lines.markersize': 3,
    'lines.markeredgewidth': 1, # capthick
    'errorbar.capsize': 2,
}


# elinewidth = 0.5
# params = {
#     'legend.fontsize': 14,
#     'figure.figsize': (10, 8),
#     'axes.labelsize': 14,
#     'axes.titlesize':14,
#     'xtick.labelsize':11,
#     'ytick.labelsize':11,
#     'lines.linewidth': 1,
#     'lines.markersize': 1,
#     'lines.markeredgewidth': 1, # capthick
#     'errorbar.capsize': 1,
# }

# params = {
#     'legend.fontsize': 24,
#     'figure.figsize': (10, 8),
#     'axes.labelsize': 24,
#     'axes.titlesize': 24,
#     'xtick.labelsize': 18,
#     'ytick.labelsize': 18,
#     'lines.linewidth': 3,
#     'lines.markersize': 12,
#     'lines.markeredgewidth': 1, # capthick
#     'errorbar.capsize': 1,
# }
plt.rcParams.update(params)

parser = argparse.ArgumentParser(description="Compare fit files with dataset.")
parser.add_argument("mfile",          help="The dataset file.")
parser.add_argument("-o", "--output", help="The output path.")
parser.add_argument("--fit", nargs="*", help="The fit files.")
args = parser.parse_args()

mfile: str = args.mfile
folder: str = os.path.dirname(mfile)

if args.output:
    folder = args.output
    if not os.path.exists(folder):
        os.mkdir(folder)

files: list[str] = args.fit
if files:
    for f in files:
        if not os.path.exists(f):
            print(f"File \"{f}\" does not exist.")
            sys.exit(1)
else:
    print("No fit files specified. Searching for fit files in the same folder as the dataset.")
    files = [f for f in os.listdir(folder) if f.endswith(('.fit', ".xvg"))]
    files = [os.path.join(folder, f) for f in files]

fits = []
labels = []
data = np.loadtxt(mfile, skiprows=1)
def load_fit(fitdata, title: str):

    # define the chi2 function
    def chi2(ymodel):
        return np.sum(((data[:, 1] - ymodel) / data[:, 2]) ** 2)

    # perform a linear fit of fitdata to the measurement
    popt, _ = curve_fit(lambda x, a, b: a*x + b, fitdata[:, 1], data[:, 1], sigma=data[:, 2], absolute_sigma=True, p0=[1, 0])
    fitdata[:, 1] = popt[0] * fitdata[:, 1] + popt[1]
    fits.append(fitdata)
    # print(f"Fit {title} to {mfile} with parameters {popt}.")

    chi2r = chi2(fits[-1][:, 1]) / (len(data[:, 1]) - 3)
    labels.append(r"$\chi^2_{red} = " + f"{chi2r:.3f}$ " + title)

# parse each file
for f in files:
    # get the stem of f without the path or extension
    stem = os.path.splitext(os.path.basename(f))[0]
    print("\tParsing file: " + stem)

    if "foxs".lower() in stem.lower():
        fitdata = np.loadtxt(f, skiprows=3, usecols=[0, 3])
        load_fit(fitdata, stem.lower())

    elif "crysol".lower() in stem.lower():
        fitdata = np.loadtxt(f, skiprows=1, usecols=[0, 3])
        load_fit(fitdata, stem.lower())

    elif "pepsi".lower() in stem.lower():
        fitdata = np.loadtxt(f, skiprows=0, comments="#", usecols=[0, 3])
        load_fit(fitdata, stem.lower())

    elif "waxsis".lower() in stem.lower():
        fitdata = np.loadtxt(f, skiprows=0, comments="#", usecols=[0, 1])
        x = fitdata[:, 0]
        y = np.interp(data[:, 0], x, fitdata[:, 1])
        fitdata = np.vstack((data[:, 0], y)).T
        load_fit(fitdata, stem.lower())

    elif "gromacs".lower() in stem.lower():
        fitdata = np.loadtxt(f, skiprows=0, comments=["@", "#", "&"], usecols=[0, 1])
        # interpolate the data to match the dataset
        x = fitdata[:, 0]/10
        y = np.interp(data[:, 0], x, fitdata[:, 1])
        fitdata = np.vstack((data[:, 0], y)).T
        load_fit(fitdata, stem.lower())

    elif "fit".lower() in stem.lower():
        fitdata = np.loadtxt(f, skiprows=1, usecols=[0, 1])
        load_fit(fitdata, stem.lower())
    else:
        print(f"Unknown fit file: \"{f}\"")

# plot the data
fig, ax = plt.subplots(figsize=(10, 6))
plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], fmt='k.', elinewidth=elinewidth, zorder=0)
for f, l in zip(fits, labels):
    plt.plot(f[:, 0], f[:, 1], label=l)

plt.xlabel("q [$Å^{-1}$]")
plt.ylabel("Intensity [arb]")
plt.title(os.path.basename(mfile.split('.')[0]))
plt.legend()
plt.semilogy()
fig.savefig(folder + '/compare_log.png', bbox_inches='tight')

plt.semilogx()
fig.savefig(folder + '/compare_loglog.png', bbox_inches='tight')

ax.set_xlim(0.1, 0.5)
fig.savefig(folder + '/compare_loglog_short.png', bbox_inches='tight')