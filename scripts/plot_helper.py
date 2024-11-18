import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from enum import Enum
import os

marker_scaling = 2

def isfloat(value: str):
    try:
        float(value)
        return True
    except ValueError:
        return False


class PlotType(Enum):
    Landscape = "PlotLandscape"
    Dataset = "PlotDataset"
    Histogram = "PlotHistogram"
    Image = "PlotImage"
    ImageAtoms = "PlotImageAtoms"
    Hline = "PlotHline"
    Vline = "PlotVline"
    Invalid = ""

# Options class. Defines how a dataset is plotted.
class Options:
    def __init__(self):
        # visuals
        self.color = "k"
        self.linestyle = "-"
        self.markerstyle = "."
        self.linewidth = 1
        self.markersize = 1
        self.drawline = True
        self.drawmarker = False
        self.drawerror = False
        self.zorder = 0

        # labels
        self.title = ""
        self.xlabel = "x"
        self.ylabel = "y"
        self.zlabel = "z"
        self.legend = ""

        # axes
        self.xrange = []
        self.yrange = []
        self.xlog = False
        self.ylog = False
        self.xshift = 0

        # other stuff
        self.dof = 0 # degrees of freedom for chi2 plots
        self.stagger = 1 # factor to multiply the y values by
        self.normalize = 0

    def parse_option(self, line):
        words = line.split()
        if len(words) < 2:
            return

        if words[0].startswith("#"):
            return

        # visuals
        match (words[0]):
            case "color":
                self.color = words[1]
            case "line_style":
                self.linestyle = words[1]
            case "marker_style":
                self.markerstyle = words[1]
            case "line_width":
                self.linewidth = float(words[1])
            case "marker_size":
                self.markersize = float(words[1])
            case "draw_line":
                self.drawline = int(words[1])
            case "draw_markers":
                self.drawmarker = int(words[1])
            case "draw_errors":
                self.drawerror = int(words[1])
            case "zorder":
                self.zorder = int(words[1])

        # labels
            case "title":
                self.title = " ".join(words[1:])
            case "xlabel":
                self.xlabel = " ".join(words[1:])
            case "ylabel":
                self.ylabel = " ".join(words[1:])
            case "zlabel":
                self.zlabel = " ".join(words[1:])
            case "legend":
                self.legend = " ".join(words[1:])
        
        # axes
            case "xlimits":
                if (words[1] == words[2]):
                    return
                self.xrange = [float(words[1]), float(words[2])]
            case "ylimits":
                if (words[1] == words[2]):
                    return
                self.yrange = [float(words[1]), float(words[2])]
            case "logx":
                self.xlog = int(words[1])
            case "logy":
               self.ylog = int(words[1])
            case "xshift":
                self.xshift = float(words[1])

        # other stuff
            case "dof":
                self.dof = int(words[1])
            case "stagger":
                self.stagger = float(words[1])
            case "normalize":
                self.normalize = float(words[1])

        # invalid option
            case _:
                print("Options.parse_option: Invalid option: " + words[0])
                exit(1)

class Dataset:
    def __init__(self, data, options):
        self.data = np.array(data)
        self.options = options

class Hline:
    def __init__(self, y, options):
        self.y = y
        self.options = options

class Vline:
    def __init__(self, x, options):
        self.x = x
        self.options = options

def read_options(file):
    """Reads plot options from a .plot file.

    Args:
        file: The file iterator.

    Returns:
        options: The plot options.
    """

    options: Options = Options()
    while(line := file.readline()):
        line = line.rstrip()

        # empty lines separates sections
        if line == "":
            break

        options.parse_option(line)
    return options

def read_dataset(file):
    """Reads a dataset from a .plot file.

    Args:
        file: The file iterator.

    Returns:
        datasets: A list of datasets.
    """

    data: list[tuple[float]] = []
    while(line := file.readline()):
        line = line.rstrip()

        # empty lines separates sections
        if line == "":
            break

        words = line.split()

        # otherwise just parse the line
        l = [float(w) for w in words]
        data.append(l)
        continue

    options = read_options(file)

    # normalize the data to start at 1
    if options.normalize != 0:
        scaling = data[0][1]/options.normalize
        if options.stagger != 1:
            scaling /= options.stagger
            options.stagger = 1

        if scaling != 0:
            for i in range(len(data)):
                data[i][1] /= scaling
            if len(data[0]) > 2:
                for i in range(len(data)):
                    data[i][2] /= scaling

    if options.stagger != 1:
        for i in range(len(data)):
            data[i][1] *= options.stagger
        if len(data[0]) > 2:
            for i in range(len(data)):
                data[i][2] *= options.stagger

    return Dataset(data, options)

def read_hline(file):
    """Reads a horizontal line from a .plot file.

    Args:
        file: The file iterator.

    Returns:
        hline: The horizontal line.
    """

    val = float(file.readline().rstrip())
    file.readline() # empty space
    return Hline(val, read_options(file))

def read_vline(file):
    """Reads a vertical line from a .plot file.

    Args:
        file: The file iterator.

    Returns:
        vline: The vertical line.
    """

    val = float(file.readline().rstrip())
    file.readline() # empty space
    return Vline(val, read_options(file))

def read_2dhist(file):
    """Reads a 2d histogram from a .plot file.

    Args:
        file: The file iterator.

    Returns:
        hist: The 2d histogram.
    """

    # read axes
    xline = file.readline().rstrip().split() # format is "xaxis: <n> <min> <max>"
    x_axis = [int(xline[1]), float(xline[2]), float(xline[3])]
    yline = file.readline().rstrip().split()
    y_axis = [int(yline[1]), float(yline[2]), float(yline[3])]

    x, y = [], []
    xstep = (x_axis[2] - x_axis[1])/x_axis[0]
    for i in range(x_axis[0]):
        x.append(x_axis[1] + i*xstep)

    ystep = (y_axis[2] - y_axis[1])/y_axis[0]
    for i in range(y_axis[0]):
        y.append(y_axis[1] + i*ystep)

    z = []
    while(line := file.readline()):
        line = line.rstrip()

        # empty lines separates sections
        if line == "":
            break

        words = line.split()
        z.append([float(w) for w in words])

    return [np.array(x), np.array(y), np.array(z), read_options(file)]

first_plot: bool = True
def plot_dataset(d: Dataset):
    """Plots a dataset.

    Args:
        data: A list of datasets to plot.
        options: Plot options.
    """

    # divide by degrees of freedom if present. Easy fix to support reduced chi2 landscapes.
    if d.options.dof != 0:
        d.data[:,1] = d.data[:,1]/d.options.dof

    if d.options.xshift != 0:
        d.data[:,0] += d.options.xshift

    if d.options.drawerror:
        if (d.data.shape[1] < 3):
            print("plot_dataset: Not enough columns for error bars.")
            exit(1)
        plt.errorbar(d.data[:,0], d.data[:,1], yerr=d.data[:,2],
            color=d.options.color,
            linestyle="none",
            marker=d.options.markerstyle,
            markersize=d.options.markersize*marker_scaling,
            label=d.options.legend,
            capsize=2*marker_scaling,
            zorder=5
        )

    if d.options.drawmarker and d.options.drawline:
        plt.plot(d.data[:,0], d.data[:,1],
            color=d.options.color,
            linestyle=d.options.linestyle,
            linewidth=d.options.linewidth*marker_scaling,
            marker=d.options.markerstyle,
            markersize=d.options.markersize*marker_scaling,
            label=d.options.legend,
            zorder=d.options.zorder
        )

    elif d.options.drawmarker:
        plt.plot(d.data[:,0], d.data[:,1],
            color=d.options.color,
            linestyle="none",
            marker=d.options.markerstyle,
            markersize=d.options.markersize*marker_scaling,
            label=d.options.legend,
            zorder=d.options.zorder
        )

    elif d.options.drawline:
        plt.plot(d.data[:,0], d.data[:,1],
            color=d.options.color,
            linestyle=d.options.linestyle,
            linewidth=d.options.linewidth*marker_scaling,
            label=d.options.legend,
            zorder=d.options.zorder
        )

    global first_plot
    if (first_plot):
        first_plot = False
        plt.title(d.options.title)
        plt.xlabel(r"{}".format(d.options.xlabel))
        plt.ylabel(r"{}".format(d.options.ylabel))
        plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 4))
        if (d.options.xrange != []):
            plt.xlim(d.options.xrange)
        if (d.options.yrange != []):
            plt.ylim(d.options.yrange)
        if (d.options.xlog):
            plt.xscale("log")
        if (d.options.ylog):
            plt.yscale("log")
    if (d.options.legend):
        plt.legend()
    return

def plot_hline(h: Hline):
    """Plots a horizontal line.

    Args:
        h: The horizontal line to plot.
    """

    plt.axhline(h.y,
        color=h.options.color,
        linestyle=h.options.linestyle,
        linewidth=h.options.linewidth*marker_scaling,
        label=h.options.legend
    )
    if (h.options.legend):
        plt.legend()
    return

def plot_vline(v: Vline):
    """Plots a vertical line.

    Args:
        v: The vertical line to plot.
    """

    plt.axvline(v.x,
        color=v.options.color,
        linestyle=v.options.linestyle,
        linewidth=v.options.linewidth*marker_scaling,
        label=v.options.legend
    )
    if (v.options.legend):
        plt.legend()
    return

def plot_landscape(d: Dataset):
    """Plots a landscape.

    Args:
        d: The dataset to plot.
    """
    
    ax = plt.gcf().add_subplot(111, projection='3d')
    x, y, z = d.data[:, 0], d.data[:, 1], d.data[:, 2]

    ax.scatter(x, y, z, c=z, cmap=mpl.colormaps["coolwarm"])
    ax.set_xlabel(r"{}".format(d.options.xlabel))
    ax.set_ylabel(r"{}".format(d.options.ylabel))
    ax.set_zlabel(r"{}".format(d.options.zlabel))
    return

def plot_image(x, y, z):
    """Plots a 2D histogram.

    Args:
        x: The x-axis data.
        y: The y-axis data.
        z: The z-axis data.
    """

    X, Y = np.meshgrid(x, y)
    plt.contourf(X, Y, z, 100, cmap=mpl.cm.coolwarm)
    plt.colorbar()
    return

def plot_file(file: str):
    """Plots a .plot file.

    Args:
        file: The input file.
    """

    print("Plotting file: " + file)
    def determine_type(line: str) -> PlotType:
        for t in PlotType:
            if line == t.value:
                return t
        return PlotType.Invalid

    global first_plot
    first_plot = True
    plt.figure()
    with open(file) as f:
        # keep reading until eof
        while line := f.readline():
            # determine the type of plot to make
            match determine_type(line.rstrip()):
                case PlotType.Dataset:
                    dataset: Dataset = read_dataset(f)
                    if dataset.data.size == 0:
                        if dataset.options.title != "":
                            title = dataset.options.title
                        elif dataset.options.legend != "":
                            title = dataset.options.legend
                        else:
                            title = "(unnammed)"
                        print(f"Skipping empty dataset \"{title}\"")
                        continue
                    plot_dataset(dataset)

                case PlotType.Hline:
                    hline: Hline = read_hline(f)
                    plot_hline(hline)

                case PlotType.Vline:
                    vline: Vline = read_vline(f)
                    plot_vline(vline)

                case PlotType.Landscape:
                    dataset: Dataset = read_dataset(f)
                    plot_landscape(dataset)

                case PlotType.Histogram:
                    dataset: Dataset = read_dataset(f)
                    plot_dataset(dataset)

                case PlotType.Image:
                    x, y, z, _ = read_2dhist(f)
                    plot_image(x, y, z)

                case PlotType.ImageAtoms:
                    print("ImageAtoms not implemented yet.")
                    exit(1)

                case _:
                    print("plot_file: Invalid plot type: " + line)
                    exit(1)

    path = file.rsplit('.', 1)[0]
    plt.tight_layout()
    plt.savefig(path, dpi=600)
    plt.close()

    # delete the .plot file
    # os.remove(file)

    return

from scipy.optimize import curve_fit
def plot_fits(ausaxs_file, fit_files, title=""):
    """
    Plots the given fit files. 
    """

    fits = []
    labels = []
    colors = []
    data = np.loadtxt(ausaxs_file, skiprows=2)
    header = open(ausaxs_file).readline().split()
    for entry in header:
        if entry.startswith("dof="):
            ausaxs_dof = len(data[:, 1]) - int(entry.split("=")[1])

    def chi2(ymodel):
        return np.sum(((data[:, 1] - ymodel) / data[:, 2]) ** 2)

    def load_fit(fitdata, title, dof = 3):
        # perform a linear fit of fitdata to the measurement
        popt, _ = curve_fit(lambda x, a, b: a*x + b, fitdata[:, 1], data[:, 1], sigma=data[:, 2], absolute_sigma=True, p0=[1, 0])
        fitdata[:, 1] = popt[0] * fitdata[:, 1] + popt[1]
        fits.append(fitdata)

        chi2r = chi2(fits[-1][:, 1]) / (len(data[:, 1]) - dof)
        labels.append(r"$\chi^2_{red} = " + f"{chi2r:.3f}$ " + title)

    # parse each file
    for f in fit_files:
        # get the stem of f without the path or extension
        stem = os.path.splitext(os.path.basename(f))[0]
        print("\tParsing file: " + stem)

        if "foxs".lower() in stem.lower():
            fitdata = np.loadtxt(f, skiprows=3, usecols=[0, 3])
            load_fit(fitdata, "FoXS")
            colors.append("tab:orange")

        elif "crysol".lower() in stem.lower():
            fitdata = np.loadtxt(f, skiprows=1, usecols=[0, 3])
            load_fit(fitdata, "CRYSOL")
            colors.append("tab:cyan")

        elif "pepsi".lower() in stem.lower():
            fitdata = np.loadtxt(f, skiprows=0, comments="#", usecols=[0, 3])
            load_fit(fitdata, "Pepsi-SAXS")
            colors.append("tab:blue")

        elif "waxsis".lower() in stem.lower():
            fitdata = np.loadtxt(f, skiprows=0, comments="#", usecols=[0, 1])
            x = fitdata[:, 0]
            y = np.interp(data[:, 0], x, fitdata[:, 1])
            fitdata = np.vstack((data[:, 0], y)).T
            load_fit(fitdata, "WAXSiS")
            colors.append("tab:purple")

        elif "waxs_final".lower() in stem.lower():
            fitdata = np.loadtxt(f, skiprows=0, comments=["@", "#", "&"], usecols=[0, 1])
            # interpolate the data to match the dataset
            x = fitdata[:, 0]/10
            y = np.interp(data[:, 0], x, fitdata[:, 1])
            fitdata = np.vstack((data[:, 0], y)).T
            load_fit(fitdata, "GROMACS")
            colors.append("tab:green")

        else:
            print(f"Unknown fit file: \"{f}\"")

    # ausaxs plot
    load_fit(data[:,[0, 3]], "AUSAXS", ausaxs_dof)
    colors.append("tab:red")

    # plot the data
    fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    plt.sca(ax[0])
    plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], fmt='k.', zorder=0)
    for f, l, c in zip(fits, labels, colors):
        plt.plot(f[:, 0], f[:, 1], label=l, color=c)
    plt.ylabel("I(q)")
    plt.legend()
    plt.semilogy()
    plt.title(title)

    plt.sca(ax[1])
    plt.axhline(0, color='k', lw=0.5)
    for f, c in zip(fits, colors):
        plt.plot(f[:, 0], (data[:, 1] - f[:, 1]) / data[:, 2], ".", color=c)
    plt.xlabel("q [$Å^{-1}$]")
    plt.ylabel("Residuals")

    plt.tight_layout()
    fig.savefig(os.path.dirname(ausaxs_file) + '/log.png', dpi=600)
    print("Plotted log.png")

    ax[0].semilogx()
    fig.savefig(os.path.dirname(ausaxs_file) + '/loglog.png', dpi=600)
    print("Plotted loglog.png")
    return
