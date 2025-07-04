# SPDX-License-Identifier: LGPL-3.0-or-later
# Author: Kristian Lytje

from multiprocessing import freeze_support
import matplotlib.pyplot as plt
import concurrent.futures
import argparse
import sys
import os

from plot_helper import *

def parse_args():
    parser = argparse.ArgumentParser(prog="plot", description="AUSAXS plotting utility.")
    parser.add_argument("--folder", type=str, default=".", help="Folder to search for .plot files.")
    parser.add_argument("--max_depth", type=int, default=4, help="Maximum depth to search for .plot files.")
    parser.add_argument("--max_files", type=int, default=30, help="Maximum number of files to process.")
    parser.add_argument("--big", action="store_true", help="Use big font size.")
    parser.add_argument("--medium", action="store_true", help="Use medium font size.")
    parser.add_argument("--title", type=str, help="Title for the plot.")
    parser.format_help
    return parser.parse_args()

def main():
    args = parse_args()
    folder = args.folder
    max_depth = args.max_depth
    max_files = args.max_files
    if args.big:
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
    elif args.medium:
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
    else:
        params = {
            'legend.fontsize': 14,
            'figure.figsize': (10, 8),
            'axes.labelsize': 14,
            'axes.titlesize': 14,
            'xtick.labelsize': 11,
            'ytick.labelsize': 11,
            'lines.markersize': 5,
            'backend': 'Agg'
        }
    title = args.title

    plt.rcParams.update(params)

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