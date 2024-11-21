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

data_ubiq = np.loadtxt("temp/debye_waller_fitting_ubiq.txt", skiprows=1)
data_lyz = np.loadtxt("temp/debye_waller_fitting_lyz.txt", skiprows=1)

plt.figure(figsize=(16, 10))
plt.title("Debye-Waller exv factor")
plt.ylabel("Count")
plt.xlabel("Value")
plt.hist(data_ubiq[:, 0], bins=25, alpha=0.5, label='ubiquitin')
plt.hist(data_lyz[:, 0], bins=25, alpha=0.5, label='lysozyme')
plt.legend()
plt.show()

plt.figure(figsize=(16, 10))
plt.title("Debye-Waller atomic factor")
plt.ylabel("Count")
plt.xlabel("Value")
plt.hist(data_ubiq[:, 1], bins=25, alpha=0.5, label='ubiquitin')
plt.hist(data_lyz[:, 1], bins=25, alpha=0.5, label='lysozyme')
plt.legend()
plt.show()

plt.figure(figsize=(16, 10))
plt.title("Debye-Waller exv factor (no atomic)")
plt.ylabel("Count")
plt.xlabel("Value")
plt.hist(data_ubiq[:, 2], bins=25, alpha=0.5, label='ubiquitin')
plt.hist(data_lyz[:, 2], bins=25, alpha=0.5, label='lysozyme')
plt.legend()
plt.show()