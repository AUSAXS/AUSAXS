import matplotlib.pyplot as plt
import numpy as np
import sys
import os

from plot_helper import *

# handle command line arguments
help = "Usage: plot <folder>"
params = {
    'legend.fontsize': 14,
    'figure.figsize': (10, 8),
    'axes.labelsize': 14,
    'axes.titlesize':14,
    'xtick.labelsize':11,
    'ytick.labelsize':11
}

match len(sys.argv):
    case 1: 
        if os.path.exists("figures"):
            folder = "figures"
        else: 
            print(help)
            exit(0)
    case 2:
        if (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
            print("Plotting tool for the .plot files from AUSAXS.")
            print(help)
            exit(0)
        else:
            folder = sys.argv[1]
            if not os.path.exists(folder):
                print(f"Folder {sys.argv[1]} does not exist.")
                exit(1)
    case 3: # secret option to change matplotlib parameters
        folder = sys.argv[1]
        if not os.path.exists(folder):
            print(f"Folder {sys.argv[1]} does not exist.")
            exit(1)
        if sys.argv[2] == "--big":
            params = {
                'legend.fontsize': 28,
                'figure.figsize': (10, 8),
                'axes.labelsize': 28,
                'axes.titlesize': 28,
                'xtick.labelsize': 20,
                'ytick.labelsize': 20
            }
            
    case _:
        print(help)
        exit(0)

plt.rcParams.update(params)

if folder == "":
    folder = sys.argv[1]
for currentpath, folders, files in os.walk(folder):
    for file in files:
        if file.endswith(".plot"):
            print("Plotting file: " + file)
            plot_file(os.path.join(currentpath, file))