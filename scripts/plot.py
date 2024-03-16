import matplotlib.pyplot as plt
import numpy as np
import sys
import os

from plot_helper import *

# handle command line arguments
help = "Usage: plot <folder>\nPlots all .plot files in the given folder and subfolders."
params = {
    'legend.fontsize': 14,
    'figure.figsize': (10, 8),
    'axes.labelsize': 14,
    'axes.titlesize':14,
    'xtick.labelsize':11,
    'ytick.labelsize':11
}

title=""
match len(sys.argv):
    case 1: 
        if os.path.exists("output"):
            folder = "output"
        else:
            folder = "."

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
                'ytick.labelsize': 20,
                'lines.markersize': 12
            }
        elif sys.argv[2] == "--medium":
            params = {
                'legend.fontsize': 24,
                'figure.figsize': (10, 8),
                'axes.labelsize': 24,
                'axes.titlesize': 24,
                'xtick.labelsize': 18,
                'ytick.labelsize': 18,
                'lines.markersize': 10
            }
    
    case 4:
        folder = sys.argv[1]
        if sys.argv[2] == "--logtitle":
            title=sys.argv[3]
            
    case _:
        print(help)
        exit(0)

plt.rcParams.update(params)

if folder == "":
    folder = sys.argv[1]

for currentpath, folders, files in os.walk(folder):
    fit_file = ""
    dat_file = ""
    report_file = ""
    for file in files:
        extension = file.split(".")[-1]
        match extension:
            case "fit":
                fit_file = os.path.join(currentpath, file)
            case "scat":
                dat_file = os.path.join(currentpath, file)
            case "plot":
                print("Plotting file: " + file)
                plot_file(os.path.join(currentpath, file))
        
        if file == "report.txt":
            report_file = os.path.join(currentpath, file)
    
    if fit_file != "" and dat_file != "" and report_file != "":
        plot_intensity_fit(dat_file, fit_file, report_file, title)