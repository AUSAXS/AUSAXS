import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from plot_helper import read_dataset

params = {
    'legend.fontsize': 14,
    'figure.figsize': (10, 8),
    'axes.labelsize': 14,
    'axes.titlesize':14,
    'xtick.labelsize':11,
    'ytick.labelsize':11,
    'lines.linewidth': 1,
    'lines.markersize': 1,
    'lines.markeredgewidth': 1, # capthick
    'errorbar.capsize': 1,
}

# params = {
#     'legend.fontsize': 28,
#     'figure.figsize': (10, 8),
#     'axes.labelsize': 28,
#     'axes.titlesize': 28,
#     'xtick.labelsize': 20,
#     'ytick.labelsize': 20,
#     'lines.linewidth': 3,
#     'lines.markersize': 12,
#     'lines.markeredgewidth': 1, # capthick
#     'errorbar.capsize': 1,
# }
plt.rcParams.update(params)

if (len(sys.argv) != 4):
    print("Usage: python plot_intensityfit.py <data> <fit> <# to subtract from dof>")
    exit(1)

data = np.loadtxt(sys.argv[1], skiprows=1)
if sys.argv[2].endswith('.dat') or sys.argv[2].endswith('.fit'):
    fit = np.loadtxt(sys.argv[2], skiprows=1, usecols=[0, 1])
elif sys.argv[2].endswith('.plot'):
    with open(sys.argv[1]) as f:
        f.readline() # skip first line
        fit = read_dataset(f).data

dof = int(sys.argv[3])

# calculate chi2
def chi2(ymodel):
    return np.sum(((data[:, 1] - ymodel) / data[:, 2]) ** 2)

chi2r = chi2(fit[:, 1]) / (len(data[:, 1]) - dof)

# plot the data in loglog and with residuals underneath

fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
ax[0].errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], fmt='k.', zorder=0)
ax[0].plot(fit[:, 0], fit[:, 1], label=r"$\chi^2_{red} = " + f"{chi2r:.3f}$", color='red')
ax[0].set_ylabel("I(q)")
ax[0].legend()
ax[0].semilogy()
ax[0].set_title(os.path.basename(sys.argv[1].split('.')[0]))

ax[1].axhline(0, color='k', lw=0.5)
ax[1].plot(data[:, 0], (data[:, 1] - fit[:, 1]) / data[:, 2], 'k.')
ax[1].set_xlabel("q")
ax[1].set_ylabel("Residuals")

fig.savefig(os.path.dirname(sys.argv[1]) + '/log.png', dpi=600)
print("Plotted log.png")

ax[0].semilogx()
fig.savefig(os.path.dirname(sys.argv[1]) + '/loglog.png', dpi=600)
print("Plotted loglog.png")