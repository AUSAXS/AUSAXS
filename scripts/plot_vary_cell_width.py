from plot_extensions import load_file 
from plot_helper import *

import matplotlib.pyplot as plt

data = load_file("output/vary_grid_cell_size/intensities.png.plot")

fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
plt.sca(ax[0])
for d in data[:4]:
    if isinstance(d, Dataset):
        plot_dataset(d)

ax[0].set_xscale("linear")
ax[0].set_ylim(1e-4, 1e-2)

plt.sca(ax[1])
data = load_file("output/vary_grid_cell_size/difference.png.plot")
for d in data[:4]:
    if isinstance(d, Dataset):
        plot_dataset(d)

ax[1].set_xscale("linear")
ax[1].set_ylim(0.5, 1.05)
plt.show()