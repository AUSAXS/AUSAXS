import matplotlib.pyplot as plt
import numpy as np
from enum import Enum
import os

def isfloat(value: str):
    try:
        float(value)
        return True
    except ValueError:
        return False


class PlotType(Enum):
    Dataset = "PlotDataset"
    Histogram = "PlotHistogram"
    Image = "PlotImage"
    ImageAtoms = "PlotImageAtoms"
    Hline = "PlotHline"
    Vline = "PlotVline"

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

        # labels
        self.title = ""
        self.xlabel = "x"
        self.ylabel = "y"
        self.legend = ""

        # axes
        self.xrange = []
        self.yrange = []
        self.xlog = False
        self.ylog = False

    def parse_option(self, line):
        words = line.split()
        if len(words) < 2:
            return

        # visuals
        if (words[0] == "color"):
            self.color = words[1]
        elif (words[0] == "line_style"):
            self.linestyle = words[1]
        elif (words[0] == "marker_style"):
            self.markerstyle = words[1]
        elif (words[0] == "line_width"):
            self.linewidth = float(words[1])
        elif (words[0] == "marker_size"):
            self.markersize = float(words[1])
        elif (words[0] == "draw_line"):
            self.drawline = int(words[1])
        elif (words[0] == "draw_markers"):
            self.drawmarker = int(words[1])
        elif (words[0] == "draw_errors"):
            self.drawerror = int(words[1])

        # labels
        elif (words[0] == "title"):
            self.title = words[1]
        elif (words[0] == "xlabel"):
            self.xlabel = words[1]
        elif (words[0] == "ylabel"):
            self.ylabel = words[1]
        elif (words[0] == "legend"):
            self.legend = words[1]
        
        # axes
        elif (words[0] == "xlimits"):
            if (words[1] == words[2]):
                return
            self.xrange = [float(words[1]), float(words[2])]
        elif (words[0] == "ylimits"):
            if (words[1] == words[2]):
                return
            self.yrange = [float(words[1]), float(words[2])]
        elif (words[0] == "logx"):
            self.xlog = int(words[1])
        elif (words[0] == "logy"):
            self.ylog = int(words[1])
        elif (words[0] == "same"):
            self.same = int(words[1])
        else:
            print("Error: Invalid option: " + words[0])
            exit(1)

class Dataset:
    def __init__(self, data, options):
        self.data = np.array(data)
        self.options = options

def read_dataset(file: str):
    """Reads a dataset from a .plot file.

    Args:
        file: The file object to read from.

    Returns:
        datasets: A list of datasets.
    """

    datasets: list[Dataset] = []
    with open(file, "r") as f:
        data: list[tuple[float]] = []
        options: Options = Options()
        section: int = 0
        f.readline() # skip first line
        for line in f:
            line = line.rstrip()

            # empty lines separates sections
            if line == "":
                section += 1
                continue

            match section:
                # SECTION 0: reading the dataset itself
                case 0:
                    words = line.split()

                    # otherwise just parse the line
                    l = [float(w) for w in words]
                    data.append(l)
                    continue

                # SECTION 1: reading the plot options
                case 1: 
                    options.parse_option(line)
                    continue

                # SECTION 2: checking for more plots
                case 2: 
                    if line == PlotType.Dataset.value:
                        datasets.append(Dataset(data, options))
                        data = []
                        options = Options()
                        section = 0
                        continue
                case _:
                    break

    datasets.append(Dataset(data, options))
    return datasets

def plot_datasets(datasets: list[Dataset], path: str): 
    """Plots a dataset.

    Args:
        data: A list of datasets to plot.
        options: Plot options.
    """

    plt.figure(figsize=(10, 8))
    first: bool = True
    for d in datasets:
        if d.options.drawerror:
            if (d.data.shape[1] < 3):
                print("Error: Not enough columns for error bars.")
                exit(1)
            plt.errorbar(d.data[:,0], d.data[:,1], yerr=d.data[:,2], 
                color=d.options.color, 
                linestyle="none",
                marker=d.options.markerstyle, 
                markersize=d.options.markersize, 
                capsize=2, 
            )

        elif d.options.drawmarker: 
            plt.plot(d.data[:,0], d.data[:,1], 
                color=d.options.color, 
                linestyle="none",
                marker=d.options.markerstyle, 
                markersize=d.options.markersize, 
            )

        if d.options.drawline:
            plt.plot(d.data[:,0], d.data[:,1], 
                color=d.options.color, 
                linestyle=d.options.linestyle, 
                linewidth=d.options.linewidth, 
            )
        
        if (first):
            first = False
            plt.title(d.options.title)
            plt.xlabel(d.options.xlabel)
            plt.ylabel(d.options.ylabel)
            if (d.options.xrange != []):
                plt.xlim(d.options.xrange)
            if (d.options.yrange != []):
                plt.ylim(d.options.yrange)
            if (d.options.xlog):
                plt.xscale("log")
            if (d.options.ylog):
                plt.yscale("log")
            if (d.options.legend != ""):
                plt.legend()
        
    path = path.rsplit('.', 1)[0]
    plt.savefig(path)
    return

def plot_file(file: str):
    """Plots a .plot file.

    Args:
        file: The input file.
    """

    def determine_type(file: str):
        """Determines the type of the file.

        Args:
            file: The input file.

        Returns:
            The PlotType of the file.
        """

        with open(file) as f:
            line = f.readline().rstrip()
            for p in PlotType:
                if line == p.value:
                    return p

    match determine_type(file):
        case PlotType.Dataset:
            datasets: list[Dataset] = read_dataset(file) # read the datasets
            path: str = file.rsplit('.', 1)[0] # remove the .plot extension
            plot_datasets(datasets, path)  # actually plot the datasets

        case PlotType.Histogram:
            print("Histogram not implemented yet.")
            exit(1)

        case PlotType.Image:
            print("Image not implemented yet.")
            exit(1)

        case PlotType.ImageAtoms:
            print("ImageAtoms not implemented yet.")
            exit(1)

        case PlotType.Hline:
            print("Hline not implemented yet.")
            exit(1)

        case PlotType.Vline:
            print("Vline not implemented yet.")
            exit(1)

    # delete the .plot file
    # os.remove(file)

    return