---
title: Plotting scripts - AUSAXS
description: Turn the .plot files written by AUSAXS into figures of SAXS fits and chi-square landscapes, and customize them.
---

# Plotting scripts

AUSAXS utilizes the custom plot format `.plot`. These are highly customizable text files, containing a simplified description of each part of a plot.

When using the command-line tools, you must manually convert these files to images. To do so, run the Python script `python scripts/plot.py` (found in the AUSAXS repository) from the folder containing the files. The required Python packages can be installed with `pip install matplotlib numpy scipy`. 

The script also takes an optional argument specifying the path. It will iteratively search through nested folders for `.plot` files. 

> [!NOTE]
> The script only recurses 3 folders deep, and skips folders with more than 30 items for safety.

When invoked on a directory containing an `ausaxs.fit` file, it will convert this to a semi-logarithmic and double-logarithmic plot of the fit. Note that the `ausaxs.fit` file is also compatible with the ATSAS plotting utility. 

## Customizing the plots
> [!TIP]
> (Experimental) To easily change most plot scalings, supply the optional `--medium` or `--big` parameters to the program invocation, such that the call becomes `./plot <directory> <size>`. Note that a dot `.` indicates the current directory.

Minor modifications can easily be made on the plots. To understand how, we first have to understand the file format itself.  

The file itself contains a number of curve data sections, each starting with an identifier. A typical curve is stored as:

```
PlotDataset
x1 y1 yerr1
x2 y2 yerr2
...

PlotOptions
color tab:orange
logx 1
logy 0
...
```
Here, the `PlotDataset` section contains the data itself, while the `PlotOptions` section contains information on how it should be plotted. To change e.g. the color of a curve, simple replace "tab:orange" with some other color [supported by matplotlib](https://matplotlib.org/stable/gallery/color/named_colors.html#base-colors). The title, axes, legends, and so on, can all be changed in this manner. 

To add or remove curves from a given plot, simply add or remove its PlotDataset and corresponding PlotOptions data sections from the .plot file. 
