import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

params = {
    'legend.fontsize': 16,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'lines.markersize': 10
}
plt.rcParams.update(params)

data = np.loadtxt('data/rigidbody_overlap_1.txt')

xx = np.linspace(0, 3, 100)
data = data[data[:,0] < xx[-1],:]
data[:,1:] /= np.amax(data[:,1:])

plt.figure(figsize=(10, 6))
plt.hist(data[:,0], bins=data[:,0], weights=data[:,1], histtype='step', label=r'current conformation', lw=2)
plt.hist(data[:,0], bins=data[:,0], weights=data[:,2], histtype='step', label=r'initial conformation', lw=2)
plt.plot(xx, np.exp(-5*xx), "k", label='$\exp(-5x)$', lw=2)
plt.xlim([xx[0], xx[-1]])
plt.ylim([7e-5, 2])
plt.semilogy()
plt.legend()
plt.ylabel("Normalized count")
plt.xlabel("Distance [Å]")
plt.tight_layout()
plt.savefig("rigidbody_overlap.png", dpi=300)
# plt.show()
