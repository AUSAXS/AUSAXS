import numpy as np
import matplotlib.pyplot as plt

params = {
    'legend.fontsize': 24,
    'figure.figsize': (10, 8),
    'axes.labelsize': 24,
    'axes.titlesize': 24,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'lines.markersize': 10
}
plt.rcParams.update(params)

data = np.loadtxt("temp/debye_waller_fit_vals.dat")

plt.figure(figsize=(16, 10))
plt.title("Debye-Waller exv factor")
plt.ylabel("Count")
plt.xlabel("Value")
plt.hist(data[:, 0], bins=25, alpha=1)
plt.show()

plt.figure(figsize=(16, 10))
plt.title("Debye-Waller atomic factor")
plt.ylabel("Count")
plt.xlabel("Value")
plt.hist(data[:, 1], bins=25, alpha=1)
plt.show()