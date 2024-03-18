import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np
import sys
import os
from plot_helper import read_dataset

img_scaling = 1
params = {
    'legend.fontsize': 14,
    'figure.figsize': (10, 8),
    'axes.labelsize': 14,
    'axes.titlesize':14,
    'xtick.labelsize':11,
    'ytick.labelsize':11,
    'lines.linewidth': 1.5,
    'lines.markersize': 1,
    'lines.markeredgewidth': 1, # capthick
    'errorbar.capsize': 2,
    'lines.markersize': 5
}

#img_scaling = 0.9
# params = {
#     'legend.fontsize': 20,
#     'figure.figsize': (10, 8),
#     'axes.labelsize': 20,
#     'axes.titlesize': 20,
#     'xtick.labelsize': 16,
#     'ytick.labelsize': 16,
#     'lines.markersize': 7.5,
#     'lines.linewidth': 1.5
# }
plt.rcParams.update(params)

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
                'lines.markersize': 10,
                'lines.linewidth': 2
            }
        elif sys.argv[2] == "--medium":
            params = {
                'legend.fontsize': 24,
                'figure.figsize': (10, 8),
                'axes.labelsize': 24,
                'axes.titlesize': 24,
                'xtick.labelsize': 18,
                'ytick.labelsize': 18,
                'lines.markersize': 7.5,
                'lines.linewidth': 1.5
            }
    
    case 4:
        folder = sys.argv[1]
        if sys.argv[2] == "--title":
            title=sys.argv[3]
            
    case _:
        print(help)
        exit(0)

def plot_intensity_fit(data_file, fit_file, report_file, pymol_file, title=""):
    """Plots an intensity fit."""
    data = np.loadtxt(data_file, skiprows=1)
    fit = np.loadtxt(fit_file, skiprows=1, usecols=[0, 1])

    # calculate dof
    dof = 0
    if report_file != "":
        # find dof in report
        with open(report_file, "r") as f:
            par_section = False
            for line in f:
                if "PAR" in line:
                    par_section = True
                    continue
                if par_section:
                    if "+----" in line:
                        break
                    dof += 1

    # calculate chi2
    def chi2(ymodel):
        return np.sum(((data[:, 1] - ymodel) / data[:, 2]) ** 2)

    print(f"CHI2 is {chi2(fit[:, 1])} for {len(data[:, 1])} data points and {dof} degrees of freedom.")
    chi2r = chi2(fit[:, 1]) / (len(data[:, 1]) - dof)

    # plot the data in loglog and with residuals underneath

    fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    ax[0].errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], fmt='k.', zorder=0)
    ax[0].plot(fit[:, 0], fit[:, 1], label=r"$\chi^2_{red} = " + f"{chi2r:.3f}$", color='red')
    ax[0].set_ylabel("I(q)")
    ax[0].legend()
    ax[0].semilogy()
    if title != "":
        ax[0].set_title(title)
    else: 
        ax[0].set_title(os.path.basename(os.path.abspath(data_file).split('.')[0]))

    ax[1].axhline(0, color='k', lw=0.5)
    ax[1].plot(data[:, 0], (data[:, 1] - fit[:, 1]) / data[:, 2], 'k.')
    ax[1].set_xlabel("q [$Å^{-1}$]")
    ax[1].set_ylabel("Residuals")

    plt.tight_layout()
    fig.savefig(os.path.dirname(data_file) + '/log.png', dpi=600)
    print("Plotted log.png")

    # insert small pymol image of the structure
    imarray = plt.imread(pymol_file)
    imagebox = OffsetImage(imarray, zoom=0.22*img_scaling)
    ab = AnnotationBbox(imagebox, (0.015, 0.03), xycoords='axes fraction', box_alignment=(0, 0), frameon=False)
#    ax[0].text(0.11, 0.03, "fitted structure", transform=ax[0].transAxes, fontsize=12)
    ax[0].add_artist(ab)

    ax[0].semilogx()
    fig.savefig(os.path.dirname(data_file) + '/loglog.png', dpi=600)
    print("Plotted loglog.png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], fmt='k.', zorder=0)
    ax.plot(fit[:, 0], fit[:, 1], label=r"$\chi^2_{red} = " + f"{chi2r:.3f}$", color='red')
    ax.set_ylabel("I(q)")
    ax.set_xlabel("q [$Å^{-1}$]")
    ax.legend()
    ax.semilogy()
    if title != "":
        ax.set_title(title)
    else: 
        ax.set_title(os.path.basename(os.path.abspath(data_file).split('.')[0]))
    plt.tight_layout()
    imagebox = OffsetImage(imarray, zoom=0.2*img_scaling)
    ab = AnnotationBbox(imagebox, (0.015, 0.03), xycoords='axes fraction', box_alignment=(0, 0), frameon=False)
#    ax.text(0.11, 0.03, "fitted structure", transform=ax.transAxes, fontsize=12)
    ax.add_artist(ab)
    ax.semilogx()
    fig.savefig(os.path.dirname(data_file) + '/loglog_no_residuals.png', dpi=600)
    print("Plotted loglog_no_residuals.png")
    return

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
        
        if file == "report.txt":
            report_file = os.path.join(currentpath, file)
    
    if fit_file != "" and dat_file != "" and report_file != "":
        plot_intensity_fit(dat_file, fit_file, report_file, os.path.join(currentpath, "pymol.png"), title)