import matplotlib.pyplot as plt
import numpy as np

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


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

def read_dataset(file):
    """Reads a dataset from a .plot file.

    Args:
        file: The file object to read from.

    Returns:
        data, lines iterated
    """

    data = []
    section = 0
    options = Options()
    for line in file:
        line = line.rstrip()

        # skip empty lines
        if line == "":
            continue

        # SECTION 0: reading the dataset itself
        if (section == 0):
            words = line.split()

            # check if we're done reading the dataset
            if (words[0] == "PlotOptions"): 
                section = 1
                continue

            # otherwise just parse the line
            l = [float(w) for w in words]
            data.append(l)
            continue

        # SECTION 1: reading the plot options
        if (section == 1):
            # check if we're done reading the options
            if (line == "END"): 
                break
            options.parse_option(line)
            continue
        break

    return np.array(data), options

def plot_dataset(data, options, path): 
    """Plots a dataset.

    Args:
        data: The dataset to plot.
        options: Plot options.
    """

    if options.drawerror:
        if (data.shape[1] < 3):
            print("Error: Not enough columns for error bars.")
            exit(1)
        plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], 
            color=options.color, 
            linestyle=options.linestyle, 
            marker=options.markerstyle, 
            linewidth=options.linewidth, 
            markersize=options.markersize, 
        )

    elif options.drawmarker: 
        plt.plot(data[:,0], data[:,1], 
            color=options.color, 
            linestyle=options.linestyle, 
            marker=options.markerstyle, 
            linewidth=options.linewidth, 
            markersize=options.markersize, 
        )

    if options.drawline:
        plt.plot(data[:,0], data[:,1], 
            color=options.color, 
            linestyle=options.linestyle, 
            linewidth=options.linewidth, 
        )
    
    plt.title(options.title)
    plt.xlabel(options.xlabel)
    plt.ylabel(options.ylabel)
    if (options.xrange != []):
        plt.xlim(options.xrange)
    if (options.yrange != []):
        plt.ylim(options.yrange)
    if (options.xlog):
        plt.xscale("log")
    if (options.ylog):
        plt.yscale("log")
    if (options.legend != ""):
        plt.legend()
    
    path = path.rsplit('.', 1)[0]
    print(path)
    plt.savefig(path)
    return