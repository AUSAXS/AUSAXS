import matplotlib.pyplot as plt
import numpy as np
import sys
import os

if (len(sys.argv) != 4):
    print("Usage: python plot_intensityfit.py <data> <fit> <# to subtract from dof>")
    exit(1)

data = np.loadtxt(sys.argv[1], skiprows=1)
fit = np.loadtxt(sys.argv[2], skiprows=1, usecols=[0, 1])
dof = int(sys.argv[3])

# calculate chi2
def chi2(ymodel):
    return np.sum(((data[:, 1] - ymodel) / data[:, 2]) ** 2)

chi2r = chi2(fit[:, 1]) / (len(data[:, 1]) - dof)

# plot the data in loglog and with residuals underneath

fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
ax[0].errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], fmt='k.', zorder=0, ms=1, capsize=1, elinewidth=0.5, capthick=0.5)
ax[0].plot(fit[:, 0], fit[:, 1], label=r"$\chi^2_{red} = " + f"{chi2r:.3f}$", lw=1, color='red')
ax[0].set_ylabel("I(q)")
ax[0].legend()
ax[0].semilogy()
ax[0].set_title(os.path.basename(sys.argv[1].split('.')[0]))

ax[1].axhline(0, color='k', lw=0.5)
ax[1].plot(data[:, 0], (data[:, 1] - fit[:, 1]) / data[:, 2], 'k.', ms=1)
ax[1].set_xlabel("q")
ax[1].set_ylabel("Residuals")

fig.savefig(os.path.dirname(sys.argv[1]) + '/log.png', dpi=300)

ax[0].semilogx()
fig.savefig(os.path.dirname(sys.argv[1]) + '/loglog.png', dpi=300)