import matplotlib.pyplot as plt
import sys
import os

from plot_helper import read_dataset

# handle command line arguments
params = {
    'legend.fontsize': 14,
    'figure.figsize': (10, 8),
    'axes.labelsize': 14,
    'axes.titlesize':14,
    'xtick.labelsize':11,
    'ytick.labelsize':11
}

if len(sys.argv) != 3:
    print("Usage: plot_solv_subtraction.py <protein> <solvent>")
    exit(1)

with open(sys.argv[1]) as f:
    f.readline() # skip first line
    protein = read_dataset(f)

with open(sys.argv[2]) as f:
    f.readline() # skip first line
    solvent = read_dataset(f)

# subtract solvent from protein
subtracted = protein.data.copy()
subtracted[:,1] = protein.data[:,1] - solvent.data[:,1]

# plot
plt.figure()
plt.plot(protein.data[:,0], protein.data[:,1], label="protein")
plt.plot(solvent.data[:,0], solvent.data[:,1], label="solvent")
plt.plot(subtracted[:,0], subtracted[:,1], label="subtracted")
plt.xlabel("$q Å^-1$")
plt.ylabel("Intensity [arb]")
plt.legend()
plt.semilogy()
plt.savefig(os.path.dirname(sys.argv[1]) + "/solv_subtraction.png")