import matplotlib.pyplot as plt
import numpy as np
import sys

from plot_helper import *

# handle command line arguments
if len(sys.argv) < 2:
    print("Usage: python plot.py <filename>")
    sys.exit(1)

if (sys.argv[1] == "-h") or (sys.argv[1] == "--help"):
    print("Plotting tool for the .plot files from AUSAXS.")
    sys.exit(1)

file = sys.argv[1]
if (file[-5:] != ".plot"):
    print("Error: File must be a .plot file.")
    sys.exit(1)

# iterate through the file
with open(file, "r") as f:
    for line in f:
        line = line.rstrip()
        # skip empty lines
        if line == "":
            continue

        print(line)
        if line == "PlotDataset":
            data, options = read_dataset(f)
            print(data)
            plot_dataset(data, options, file)