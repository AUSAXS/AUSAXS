---
title: Cryo-EM map validation options - AUSAXS
description: "ausaxs em: validate cryo-EM maps (MRC/CCP4) against solution SAXS data. Covers the chi-square landscapes, output files and all arguments."
---

# Validating EM maps with SAXS data

This page documents the EM validation method: its output and all of the arguments and settings available for the fitting procedure. For a step-by-step walkthrough on real data, see the [EM validation tutorial](https://sastutorials.org/SAXS_EM/SAXS_EM.html). 

All AUSAXS tools are run through the single `ausaxs` executable, with the name of the tool as the first argument. The EM tool requires two arguments: a location of a map file, and a location of a SAXS data file:  
`ausaxs em <map file> <saxs file>`

If you leave out the SAXS data file and instead give a threshold level with `em --level <σ>`, the program runs in _simulation mode_: the scattering profile of the map at that threshold is calculated and saved as `<map-name>.dat`, without any fitting. 

## Output
By default all output will be saved in a new `output/em/<data-name>/<map-name>/` folder, although this location can be changed by specifying the `--output` argument (note that the subfolders will always be created). In simulation mode, the output is written to `output/em/simulated/<map-name>/` instead. 

In this folder you will find several `.plot` files. Turn them into images by running `ausaxs plot` in this folder. See the [CLI plotting](plotting-scripts.md) page for more information. 

After you have plotted the files, you will have the three following $\chi^2$-landscape figures:

- **chi2_evaluated_points_full**: This figure shows the full $\chi^2$-landscape for the entire range of scanned threshold values. 
- **chi2_evaluated_points_limited**: This figure shows the same range as the first, but with very high- $\chi^2$ points removed. This is done so the more interesting low- $\chi^2$ range can better be visualized. A rolling average is performed, as shown by the solid red line. This average is also the basis for checking for the presence of more local minima. All minima are shown with a dashed vertical red line. You can find dummy models for each of these minima in the `models` folder.
- **chi2_near_minimum**: This figure shows the area near the absolute minimum, sampled with the smallest step-size possible for your map. Shown are also the mean and standard deviation of the entire range. This mean is the reported average $\chi^2$ value, $\bar{\chi}^2$. The blue dot is just the interpolated $\chi^2$ value from earlier, and may thus often be somewhat higher/lower in value than the rest of the plot. 

The second and third of these plots also have an accompanying `_mass` version, where the threshold x-axis has been replaced with a mass-axis. 

Beyond these landscape plots, you will also find two figures:

- **log**: The fit and its residuals, with a linear x-axis. 
- **loglog**: Same as before, but with a logarithmic x-axis. 

The `ausaxs.fit` file contains both the data and the fitted scattering profile. The other text file is the `report.txt`, which contains information about the fit itself. Here you will find the actual parameter values, along with the fitted threshold value and the $\chi^2$. 

> [!TIP]
> The `.fit` file is compatible with the plotting utility from the ATSAS package. 

Finally you will also find the `model.pdb` file, which contains the dummy structure. Note that if your map is very big, this file may be split into several numerated parts. 

## Graphical user interface
If you prefer a visual interface, you can use the Python GUI from the [pyAUSAXS](https://github.com/AUSAXS/pyAUSAXS) package:
```bash
pip install "pyausaxs[all]"
ausaxs gui
```
If you do use this, you can skip over the next section on the command-line arguments, though not every option is available in the graphical interface.

## Arguments
Arguments are handled by [CLI11](https://github.com/CLIUtils/CLI11). You can see a list of available arguments by providing the `--help` argument: `ausaxs em --help`. 

As for the [SAXS fitter](fit.md#argument-groups), most arguments are sorted into argument _groups_. To set an argument from a group, write the group name followed by the argument, e.g. `ausaxs em <map> <data> em --levelmin 2 --frequency 2 data --qmax 0.4`. 

### Base group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`-h`, `--help`      | Help                       | Print a help message |
|`--licence`         | Print the licence          | Print the licence    |
|`-v`, `--version`   | See the program version    | See the program version. | 
|`-o`, `--output`    | Set output directory       | Set the path to the output directory. |
|`-t`, `--threads`   | Set the number of threads  | Use this to limit the number of threads. |
|`--gpu`, `--no-gpu` | Toggle GPU acceleration    | Use GPU acceleration if available. On by default once a GPU is chosen with [`ausaxs gpu install`](../../getting-started/gpu.md), unless turned off with `ausaxs gpu disable`. |

### `config` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`-f`, `--file`      | Set settings file          | Read a settings file. This file can be used to easily set new default values for most arguments. |
|`--save`            | Save settings              | Write the current settings to `settings.txt` in the current directory. |
|`--log`             | Enable logging             | Write a log file. |

### `data` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
| `--qmin`           | Set minimum q              | q-points in scattering data less than this value will be ignored. |
| `--qmax`           | Set maximum q              | q-points in scattering data larger than this value will be ignored. |
| `-u`, `--unit`     | Set the data q unit        | Change the unit of the data, either `A` or `nm`. By default an educated guess will be made. | 
| `--skip`           | Skip `N` first data points | The first `N` points in the scattering data will be skipped. |
| `--rebin`          | Rebin the data             | Rebin the data to improve the information content of each data point. |

### `em` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
| `--levelmin`, `--level` | Set minimum RMS-level | This determines the minimum threshold level used in the fit. In simulation mode, this is the threshold level the profile is calculated at. |
| `--levelmax`       | Set maximum RMS-level      | This determines the maximum threshold level used in the fit. |
| `--charge-levels`  | Set the number of charge levels | Set the number of partial histograms to use. These will be uniformly distributed between the `levelmin` and `levelmax` levels. This is primarily useful if the scanned range is very large. Note: increasing this may significantly increase the RAM usage. |
| `--frequency`      | Set the sampling frequency | This determines how often the map will be sampled for dummy atoms. |
| `--hydrate`, `--no-hydrate` | Toggle hydration  | Decide whether the dummy models will be hydrated. |
| `--fixed-weight`, `--dynamic-weight` | Toggle dynamic weights | With `--dynamic-weight`, the weight of all dummy atoms will be extracted from the map. |

### `fit` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
| `--max-iterations` | Set the maximum number of evaluations | This is a soft limit on the number of evaluations made during the fit. Typically up to 50% more evaluations than this value are used. |
| `--verbose`, `--quiet` | Toggle progress information | With `--quiet`, no optional progress information will be printed to the console. |

### Settings file
The settings file is a handy tool for setting new default values. All arguments (and even some additional options) can be set in this file, and then automatically be used in future runs. The settings file can either be manually set using the `config --file` argument, or by placing a file named `settings.txt` in the same folder as the SAXS datafile, in which case it will be automatically found and used. Options from the command line take precedence over the settings file. To get a template, run `ausaxs em config --save`. 

The most useful option beyond the arguments from the argument list is `format`, which sets the format of the outputted figures. 

### q-range
The two arguments `--qmin` and `--qmax` specify the range of q-values used for the fit. Any data points outside this range in the scattering files will be discarded. The default range is $[10^{-4},\ 0.5]\ Å^{-1}$. 

> [!IMPORTANT]
> The lower and upper limits are $[10^{-4},\ 1]\ Å^{-1}$. Values outside this range are currently not supported. 

### level-range
The two arguments `em --levelmin` and `em --levelmax` specify the range of threshold values in terms of the RMS. These are the same levels as those used by e.g. PyMOL, meaning you can use PyMOL for estimating the bounds, and then apply them with these two arguments. Note that the evaluation of higher threshold values is significantly faster than smaller values, meaning the choice of `levelmin` has a much higher impact on the performance than `--levelmax`.

> [!WARNING]
> Using too small threshold values may result in the calculations never finishing. Only lower `--levelmin` if you are absolutely certain that you need to.

### Map sampling frequency
The argument `em --frequency` determines how often dummy atoms will be sampled from the grid intrinsic in the map. A value of `2` means that only every second grid point along each axis is used, effectively reducing the size of the dummy model by a factor 8. For somewhat larger high-resolution maps, using `--frequency 2` may be enough for a decent quick fit. 

Note that since this sampling frequency directly relates to the information content in the dummy model, using higher values introduces more noise to the $\chi^2$-landscape. 

> [!TIP]
> This option is very useful for reducing computation times when combined with low threshold values where the decreased information content does not really matter. 

### To hydrate or not to hydrate
By default all dummy structures are hydrated. While our hydration algorithm runs in linear time, it will have a small impact on the runtime. If you believe your structure does not need a hydration layer, this feature can be disabled by adding the `em --no-hydrate` argument. Note that we have found the addition of this shell to essentially always be necessary, so you should only disable this feature if you are absolutely sure you do not need it. 

### Using dynamic weights
By default all dummy atoms are weighted with the same value of unity. If you would instead like to use the density information from the map itself as weights instead, you can supply the `em --dynamic-weight` argument. Note that with dynamic weights low-density noise does not contribute significantly to the scattering profile, and will thus sometimes be included in the optimized structure. See the article for more details on this. 
