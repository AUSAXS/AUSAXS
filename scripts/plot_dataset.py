import matplotlib.pyplot as plt
import numpy as np
import sys
import os

mfile = sys.argv[1]
dir = os.path.dirname(mfile)

# get all .fit files in the directory
files = [f for f in os.listdir(dir) if f.endswith('.fit')]

# load the data
data = np.loadtxt(mfile, skiprows=1)

# calculate chi2
def chi2(ymodel):
    return np.sum(((data[:, 1] - ymodel) / data[:, 2]) ** 2)

fits = []
labels = []
for f in files:
    # remove the .fit extension from f
    f = f[:-4]

    if "fit" in f:
        fits.append(np.loadtxt(dir + '/' + f + ".fit", skiprows=1, usecols=[0, 1]))
        chi2r = chi2(fits[-1][:, 1]) / (len(data[:, 1]) - 3)
        labels.append(r"$\chi^2_{red} = " + f"{chi2r:.3f}$, Fit")
    elif "foxs" in f:
        fits.append(np.loadtxt(dir + '/' + f + ".fit", skiprows=3, usecols=[0, 3]))
        chi2r = chi2(fits[-1][:, 1]) / (len(data[:, 1]) - 2)
        labels.append(r"$\chi^2_{red} = " + f"{chi2r:.3f}$, FoXS")
    elif "crysol" in f:
        fits.append(np.loadtxt(dir + '/' + f + ".fit", skiprows=1, usecols=[0, 3]))
        chi2r = chi2(fits[-1][:, 1]) / (len(data[:, 1]) - 2)
        labels.append(r"$\chi^2_{red} = " + f"{chi2r:.3f}$, Crysol")
    elif "waxsis" in f:
        fits.append(np.loadtxt(dir + '/' + f + ".fit", skiprows=0, usecols=[0, 1]))
        labels.append("WAXSiS")
    else: 
        print(f"Unknown fit file: \"{f}\"")

# plot the data
fig, ax = plt.subplots()
ax.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], fmt='k.', zorder=0, ms=1, capsize=1, elinewidth=0.5, capthick=0.5)
for f, l in zip(fits, labels):
    if ("waxsis" in str.lower(l)): continue
    ax.plot(f[:, 0], f[:, 1], label=l, lw=1)

ax.set_xlabel("q")
ax.set_ylabel("I(q)")
ax.set_title(os.path.basename(mfile.split('.')[0]))
ax.legend()
ax.semilogy()
fig.savefig(dir + '/log.pdf')

ax.semilogx()
fig.savefig(dir + '/loglog.pdf')

ax.set_xlim(0.1, 0.5)
fig.savefig(dir + '/loglog_short.pdf')