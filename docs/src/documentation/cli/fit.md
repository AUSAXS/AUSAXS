---
title: Fitting atomic models to SAXS data - AUSAXS
description: Calculate the SAXS profile of a PDB or mmCIF structure and fit it to experimental small-angle X-ray scattering data with ausaxs fit. Covers the hydration shell, excluded volume, output files and all arguments.
---

# Fitting atomic models to SAXS data

With this tool you can fit your high-resolution models to measured SAXS data. Here we will cover all of the arguments and settings which are available for the fitting procedure. 

All AUSAXS tools are run through the single `ausaxs` executable, with the name of the tool as the first argument. The most basic invocation of the fitter requires just two arguments: a path to an atomic model file, and the path to a SAXS data file:  
`ausaxs fit <model file> <saxs file>`

The order of the two files does not matter. If you leave out the SAXS data file, the program instead runs in _simulation mode_, where the scattering profile of the model is calculated without fitting it to anything. 

> [!TIP]
> `ausaxs saxs` is an alias of `ausaxs fit`.

[`ausaxs sans`](sans.md) and [`ausaxs debye`](debye.md) are built on this tool, and take most of the arguments described here. Their own pages only list where they differ.

## Output
By default all output will be saved in a new `output/saxs/<data-name>/` folder, where `<data-name>` is the name of the SAXS data file. The location can be changed with the `--output` argument (the subfolder `<data-name>` will always be created). In simulation mode, the output is written to `output/saxs/simulated/<model-name>/` instead. 

In this folder you will find the following files: 
- **ausaxs.fit**: A data file containing both the data and the fitted curve evaluated at the same `q`-values as the input SAXS data file. Note that the data may be truncated and scaled to inverse Ångström, depending on your settings. This file is compatible with the plotting utility from the ATSAS package. 
- **report.txt**: A short report containing the most important information from the fit. Here you will find: 
	- *Fevals*: The number of function evaluations.
	- *chi2*: The absolute goodness-of-fit $\chi^2$ value. 
	- *dof*: The number of degrees of freedom. Say N is the number of data points in the input SAXS file. Then this will be equal to N - P, where P is the number of fit parameters. 
	- *chi2/dof*: The reduced $\chi^2$ value, $\chi^2_r$. This is what is generally reported by other fitting programs. 
	- *fitted parameters*: A list of the fit parameters showing their optimal values & associated uncertainties. 
- **model.pdb**: The model as it was used in the fit, including the generated hydration shell. 
- **p(r).png.plot** and **profiles.png.plot**: Plot files for the distance distribution and the scattering profiles of the individual contributions. 
- **rebinned.dat**, **grid.pdb**, **exv.pdb**: Only written if requested by `data --rebin`, `grid --save`, and `exv --save` respectively. 

Running `ausaxs plot` in this folder will convert the `.plot` files to images, and also generate two figures of the fit: 
- **log.png**: This top panel of this figure shows the fit plotted on top of the input data set, with the $\chi^2_r$ shown in the legend. The bottom panel shows the residual associated with each point, which ideally should have a Gaussian distribution with a mean of 0. 
- **loglog.png**: This shows the same plots as *log.png*, except with double-logarithmic axes. 

See the [CLI plotting](plotting-scripts.md) page for more information on this step. 

## Graphical user interface

If you prefer a visual interface, you can use the Python GUI from the [pyAUSAXS](https://github.com/AUSAXS/pyAUSAXS) package:
```
pip install pyausaxs[all]
ausaxs gui
```
If you do use this, you can skip over the next section on the command-line arguments.

> [!NOTE]
> Not all options are available through this interface.

## Arguments
Arguments are handled by [CLI11](https://github.com/CLIUtils/CLI11). You can see a list of available arguments by providing the `--help` argument: `ausaxs fit --help`.  

> [!WARNING]
> Arguments must be provided in a quite atypical fashion through argument _groups_! Make sure to carefully read the following section to understand how this works.

### Argument groups
To avoid flooding the terminal with information on the many arguments the program supports, we have separated most of them into argument _groups_ (subcommands), which are all listed at the bottom of the `help` output. To get more information on the arguments contained within each, simply invoke `--help` again along with the name of a group, e.g. `ausaxs fit <group> --help`. 

To set an argument from one of these groups, simply write the group name followed by the argument. Consecutive arguments from the same group can be given without repeating the group name. 

As an example, `ausaxs fit <pdb> <data> exv --model grid --fit solv --keep data --qmax 0.4` sets the `model` argument and `fit` flag from the `exv` group, the `keep` flag from the `solv` group, and sets the `qmax` argument from the `data` group. 

### Base group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`-h`, `--help`      | Help                       | Print a help message |
|`--licence`         | Print the licence          | Print the licence    |
|`-v`, `--version`   | See the program version    | See the program version. | 
|`-o`, `--output`    | Set output directory       | Set the path to the output directory. Defaults to `output/saxs/`. |
|`--allow-unknown-atoms` | Do not exit when an unknown atom is encountered | Usually the program terminates if an unknown atom is found since its charge cannot be correctly determined. If you would like to proceed anyway, supply this argument. |
|`--allow-unknown-residues` | Do not exit when an unknown residue is encountered | As above, but for residues. |
|`--offline`         | Run in offline mode        | Prevent all network requests. Normally, definitions of unknown residues may be downloaded from the PDB. |
|`-t`, `--threads`   | Set the number of threads  | Use this to limit the number of processor threads. By default all but one thread available on the system will be used. |
|`--gpu`             | Use GPU acceleration       | Use GPU acceleration if available. |

### `config` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`-f`, `--file`      | Set settings file          | Read a settings file. This file can be used to easily set new default values for most arguments. Options given on the command line take precedence over the file. |
|`--save`            | Save settings              | Write the current settings to `settings.txt` in the current directory. If no input files are given, the program exits after writing it. |
|`--log`             | Enable logging             | Write a log file. |

If no settings file is given, the program will look for one named `settings.txt`, `setting.txt`, `setup.txt`, or `config.txt` in the same folder as the SAXS data file, and use it automatically if found. 

> [!TIP]
> To create a settings file, run e.g. `ausaxs fit exv --model fraser config --save` and then edit the resulting `settings.txt`. 

### `data` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`--qmin`           | Set minimum q              | q-points in scattering data less than this value will be ignored. |
|`--qmax`           | Set maximum q              | q-points in scattering data larger than this value will be ignored. Defaults to `0.5`. |
|`-u`, `--unit`     | Set the data unit          | Set the unit of the scattering data, either `A` or `nm`. By default an educated guess will be made. |
|`--skip`           | Skip `N` first data points | The first `N` points in the scattering data will be skipped. |
|`--rebin`          | Rebin the data             | Rebin the data to improve the information content of each data point. The rebinned data is saved as `rebinned.dat`. |
|`--weighted-bins`  | Use weighted bins          | Use weighted bins for the distance histograms. By default, this is automatically decided for each molecule. |

The `q`-range must lie within $[10^{-4},\ 1]\ Å^{-1}$. 

### `molecule` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`--use-occupancy`, `--ignore-occupancy` | Use occupancy information  | Decide whether the occupancy field in the `PDB` file is used to scale the scattering contribution of each atom. It is used by default. |

### `exv` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`-m`, `--model`     | Set excluded volume model  | Change the excluded volume model. More information on these can be found at the end of this table. |
|`--fit`             | Fit excluded volume        | Enable excluded volume fitting. This will add an additional fitting parameter for most excluded volume models. |
|`--fit-density`     | Fit solvent density        | Enable excluded solvent density fitting. This will always add an additional fitting parameter. |
|`--table`           | Set excluded volume table  | Change the displaced volume table backing the model. Options: `traube`, `vdw`, `voronoi`, `mf` (minimum fluctuation). |
|`--save`            | Save excluded volume       | Save the excluded volume as `exv.pdb`. This will only happen when using a grid-based model. |
|`--ref`, `--reference` | Excluded volume reference | Path to a reference file defining the excluded volume. |

The excluded volume model is responsible for both calculating and Fourier transforming the distance histograms. There are multiple options available for this setting.

* `simple`: The default 'simple' model, which treats all atoms as being completely equivalent. This is usually a good approximation for performing quick fits. 

* `fraser`: The Gaussian sphere-based model, as proposed by [Fraser _et al._](https://doi.org/10.1107/S0021889878014296) and first implemented in `CRYSOL`. This model uses both explicit atomic and excluded volume form factors, with the latter occupying the same positions as the real atoms. This approach is equivalent to that of `CRYSOL`, `Pepsi-SAXS`, and `FoXS`. 

* `grid`: The grid-based model, where the excluded volume is modeled as a continuous grid. This also uses explicit atomic form factors, though the grid cells share a single excluded volume form factor.  

For more information regarding how each method is fitted, see the main `AUSAXS` article. 

<details markdown="1">
<summary>Additional models</summary>

These models are not as well-tested as those listed above, and are only listed for completeness. 

* `grid-base`: An alternative grid-based model, where the excluded volume fitting simply rescales the base excluded volume curve.  

* `grid-scalable`: Another grid-based model, where the excluded volume fitting is done by "stretching" the grid. A scaling factor is therefore multiplied onto all distances between dummy exv atoms. This is very slow compared to the other grid methods.  

* `average`: This is an extension of the `simple` model, which uses individual form factors for the atomic species.  

* `crysol`, `foxs`, `pepsi`, `waxsis`: These excluded volume models implements the excluded volume as defined in their respective papers. None of them exactly reproduces the output from their respective programs.  

* `none`: Disable the excluded volume entirely. 
</details>

### `solv` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`-m`, `--model`     | Set hydration shell model  | Change the model used to generate the hydration shell. Options: `radial`, `axes`, `none`. Currently this should not be changed. |
|`--keep`, `--discard` | Keep hydration molecules   | Keep the hydration molecules present in the structure file. This will cause no further hydration molecules to be generated. By default they are discarded, and a new hydration shell is generated. |
|`--fit`, `--no-fit` | Toggle shell fitting      | Decide whether the scattering power of the hydration shell is fitted. This should currently always be enabled, as the default hydration model requires it. |

### `hydrogens` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`--keep`, `--discard` | Keep explicit hydrogens    | Keep the explicit hydrogens present in the structure file instead of accounting for them implicitly. This will also switch the excluded volume tables to a version containing non-zero hydrogen volumes. By default they are discarded. | 
|`--force-field-names` | Read force-field residue names | Read force-field residue names, such as the protonation states `HID`, `HSE` and `CYX`, or `SOL` for water, by their force-field meaning instead of their usual definitions in the PDB's Chemical Component Dictionary. |

### `grid` group
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`--rvol`            | Set min exv radius         | The volume occupied by each atom is calculated as `max(rvdw, rvol)`. This setting modifies the latter, allowing for generating either small and structured or large and homogeneous excluded volumes. |
|`-w`, `--width`     | Set the grid width         | Set the grid cell width. This will have a major impact on everything relying on the grid, including the hydration model, the volume calculations, and the grid-based excluded volume, whose dummy atoms have the same width. |
|`--save`            | Save the grid              | Save the grid as `grid.pdb`. |
