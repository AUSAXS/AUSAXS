# CLI plotting

The command-line tools write their results as text files: the fit itself as `ausaxs.fit`, and every figure as a `.plot` file, a simple description of each curve in the plot. The `ausaxs plot` tool turns these into images. The [graphical interface](fit.md#graphical-user-interface) and the [Python API](https://github.com/AUSAXS/pyAUSAXS) draw their own figures, so this is only needed for results from the command line.

`ausaxs plot` is part of the Python package, and needs its plotting extension:
```bash
pip install "pyausaxs[plots]"
```

> [!TIP]
> If you use the [precompiled executables](https://ausaxs.github.io/getting-started/precompiled-binaries/) or [compiled AUSAXS yourself](https://ausaxs.github.io/getting-started/compiling/) without the Python package, the same plotting is available as the script `scripts/plot.py` in the AUSAXS repository. It needs `pip install matplotlib numpy scipy`, and takes the same options as `ausaxs plot`, except those for animations.

## Making the figures
Run it from the output folder of a command-line tool:
```bash
ausaxs plot
```
or point it at a folder with `--folder <path>`. It searches the folder and its subfolders for results, and writes each figure as a `.png` next to the file it came from:

- Every `.plot` file becomes an image of the same name.
- An `ausaxs.fit` file becomes `log.png` and `loglog.png`, the fit on semi-logarithmic and double-logarithmic axes. Any other `.fit` or `.xvg` files in the same folder, e.g. fits from other SAXS programs, are drawn on the same axes for comparison. The `ausaxs.fit` file can also be read by the ATSAS plotting utilities.
- A numbered sequence of frames becomes an animated `.gif`.

The search stops 4 folders deep, and skips folders with more than 30 files. Change these limits with `--max_depth` and `--max_files`.

## Options
| Option | Effect |
|---|---|
| `--folder <path>` | Folder to search for results. Default: the current folder. |
| `--max_depth <n>` | How many subfolders deep to search. Default: 4. |
| `--max_files <n>` | Skip folders with more than this many files. Default: 30. |
| `--medium`, `--big` | Larger fonts, lines and markers, e.g. for slides and posters. |
| `--title <text>` | Title of the fit figures. |
| `--fps <n>` or `--duration <s>` | Frame rate, or total length, of animations. |
| `--loop <n>` | How many times animations repeat. Default: 0, forever. |

## Customizing the plots
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
Here, the `PlotDataset` section contains the data itself, while the `PlotOptions` section contains information on how it should be plotted. To change e.g. the color of a curve, simply replace "tab:orange" with some other color [supported by matplotlib](https://matplotlib.org/stable/gallery/color/named_colors.html#base-colors). The title, axes, legends, and so on, can all be changed in this manner. 

To add or remove curves from a given plot, simply add or remove its PlotDataset and corresponding PlotOptions data sections from the .plot file. 
