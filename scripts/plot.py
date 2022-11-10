import matplotlib.pyplot as plt
import numpy as np
import sys
import os

from plot_helper import *

# handle command line arguments
match len(sys.argv):
    case 1: 
        if os.path.exists("figures"):
            folder = "figures"
        else: 
            print("Usage: plot <folder>")
            exit(0)
    case 2:
        if (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
            print("Plotting tool for the .plot files from AUSAXS.")
            print("Usage: plot <folder>")
            exit(0)
        else:
            folder = sys.argv[1]
            if not os.path.exists(folder):
                print(f"Folder {sys.argv[1]} does not exist.")
                exit(1)
    case _:
        print("Usage: plot <folder>")
        exit(0)

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)

if folder == "":
    folder = sys.argv[1]
for currentpath, folders, files in os.walk(folder):
    for file in files:
        if file.endswith(".plot"):
            print("Plotting file: " + file)
            plot_file(os.path.join(currentpath, file))