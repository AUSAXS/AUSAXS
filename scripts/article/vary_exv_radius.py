import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import numpy as np

data = np.loadtxt("output/vary_exv_radius/chi2.txt", dtype=str)
r_vals = data[0, :][1:]             # extract radii
data = data[1:, :]                  # remove header
labels = data[:, 0]                 # extract labels
data = data[:, 1:].astype(float)    # remove labels and convert rest to float

# plot exv radii lines with a cmap
cmap = plt.get_cmap("autumn")
fig, ax = plt.subplots(figsize=(20, 12))
for i in range(data.shape[1]):
    if i == 0 or i == data.shape[1] - 2:
        ax.plot(labels, data[:, i], color=cmap(i / data.shape[1]), label=f"r = {r_vals[i]}")
    elif i == data.shape[1] - 1:
        ax.plot(labels, data[:, i], color="k", label=f"default")
    else:
        ax.plot(labels, data[:, i], color=cmap(i / data.shape[1]))

ax.set_xlabel("exv radius")
ax.set_ylabel("Chi2")
ax.legend()
plt.savefig("output/vary_exv_radius/chi2_lines.png", dpi=300, bbox_inches="tight")

# make imshow comparison plot
fig, ax = plt.subplots(figsize=(20, 12))
cmap = plt.get_cmap("RdYlGn_r")
bounds = [-100, -50, -40, -30, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 30, 40, 50, 100]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# create a new dataset with the difference between the default and the other radii
data_diff = np.zeros_like(data)
for i in range(data.shape[1]):
    # relative improvement
    data_diff[:, i] = 100 * (data[:, i] - data[:, -1]) / data[:, -1]
data_diff = data_diff[:, :-1]

im = ax.imshow(data_diff, interpolation='nearest', cmap=cmap, norm=norm)
for i in range(len(r_vals)-1):
    for j in range(len(labels)):
        plt.axvline(i+0.5, color="black", linewidth=1)
for i in range(len(labels)):
    plt.axhline(i+0.5, color="black", linewidth=1)

plt.xticks(np.arange(len(r_vals)-1), r_vals[:-1])
plt.yticks(np.arange(len(labels)), labels, ha="right", rotation_mode="anchor")

plt.colorbar(im, label="Relative improvement [%]", shrink=0.6)
plt.savefig("output/vary_exv_radius/chi2_imshow.png", dpi=300, bbox_inches="tight")
plt.show()