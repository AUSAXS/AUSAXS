from multiprocessing import freeze_support
import matplotlib.pyplot as plt
import concurrent.futures
import sys
import os

from plot_helper import *

max_depth = 4
max_files = 30

def main():
    help_msg = "Usage: plot <folder>\nPlots all .plot files in the given folder and subfolders."
    params = {
        'legend.fontsize': 14,
        'figure.figsize': (10, 8),
        'axes.labelsize': 14,
        'axes.titlesize':14,
        'xtick.labelsize':11,
        'ytick.labelsize':11,
        'lines.markersize': 5,
        'backend': 'Agg'
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
                print(help_msg)
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
                    'lines.markersize': 12,
                    'backend': 'Agg'
                }
            elif sys.argv[2] == "--medium":
                params = {
                    'legend.fontsize': 24,
                    'figure.figsize': (10, 8),
                    'axes.labelsize': 24,
                    'axes.titlesize': 24,
                    'xtick.labelsize': 18,
                    'ytick.labelsize': 18,
                    'lines.markersize': 10,
                    'backend': 'Agg'
                }

        case 4:
            folder = sys.argv[1]
            if sys.argv[2] == "--logtitle":
                title=sys.argv[3]
                
        case _:
            print(help_msg)
            exit(0)

    plt.rcParams.update(params)

    if folder == "":
        folder = sys.argv[1]

    def get_depth(path):
        return path.count(os.sep)

    with concurrent.futures.ProcessPoolExecutor(8) as executor:
        futures = []
        invoke_depth = get_depth(folder)
        for currentpath, _, files in os.walk(folder):
            if max_depth < get_depth(currentpath) - invoke_depth:
                continue
            if max_files < len(files):
                print(f"Skipping {currentpath} because it has too many files.")
                continue

            fit_files = []
            ausaxs_file = ""
            for file in files:
                extension = file.split(".")[-1]
                if file == "ausaxs.fit":
                    ausaxs_file = os.path.join(currentpath, file)
                    continue
                match extension:
                    case "fit" | "xvg":
                        fit_files.append(os.path.join(currentpath, file))
                    case "plot":
                        futures.append(executor.submit(plot_file, os.path.join(currentpath, file)))

            if ausaxs_file:
                futures.append(executor.submit(plot_fits, ausaxs_file, fit_files, title))

        concurrent.futures.wait(futures)

if __name__ == "__main__":
    freeze_support()
    main()