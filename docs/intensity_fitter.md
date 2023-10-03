# High-resolution SAXS fitter
With this tool you can fit your high-resolution models to measured SAXS data. Here we will cover all of the arguments and settings which are available for the fitting procedure. 
The program requires just two mandatory arguments: a path to an atomic model file, and the path to a SAXS data file:
`saxs_fitter <model file> <saxs file>`
## Output
By default all output will be saved in a new `output/saxs_fitter/<model-name>/` folder, although this location can be changed by specifying the `--output` argument (note that the subfolder `<model-name>` will always be created). 
In this folder you will find three files: 
- **\<model-name\>.dat**: A copy of the original SAXS data file, potentially scaled to different units. 
- **fit.fit**: A data file containing the fitted curve evaluated at the same `q`-values as the input SAXS data file. 
- **report.txt**: A short report containing the most important information from the fit. Here you will find: 
	- *Fevals*: The number of function evaluations.
	- *chi2*: The absolute goodness-of-fit $\chi^2$ value. 
	- *dof*: The number of degrees of freedom. Say N is the number of data points in the input SAXS file. Then this will be equal to N - P, where P is the number of fit parameters. 
	- *chi2/dof*: The reduced $\chi^2$ value, $\chi^2_r$. This is what is generally reported by other fitting programs. 
	- *fitted parameters*: You will also find a list of the fit parameters showing their optimal values & associated uncertainties. Not all fitting routines support uncertainties, so depending on the settings this may be set to 0. 
Though there are no `.plot` files, running the `plot.py` python script will still generate two figures: 
- **log.png**: This top panel of this figure shows the fit plotted on top of the input data set, with the $\chi^2_r$ shown in the legend. The bottom panel shows the residual associated with each point, which ideally should have a Gaussian distribution with a mean of 0. 
- **loglog.png**: This shows the same plots as *log.png*, except with double-logarithmic axes. 

## Arguments
Arguments are handled by [CLI11](https://github.com/CLIUtils/CLI11). You can see a list of available arguments by providing the `--help` argument. 
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`-h`, `--help`      | Help                       | Print a help message |
|`-s`, `--settings`  | Set settings file          | Set a settings file. This file can be used to easily set new default values for most agruments. |
| `-r`, `--reduce`   | Hydration reduction factor | Reduce the number of generated hydration molecules to some percentage of the number of atoms. Default: 0.1. |
|`-o`, `--output`    | Set output directory       | Set the path to the output directory |
|`-t`, `--threads`   | Set the number of threads  | Use this to limit the number of threads. By default all but one thread available on the system will be used. |
|`--quiet`           | Disbale optional information | No optional progress information will be printed to the console. |
| `--qmin`           | Set minimum q              | q-points in scattering data less than this value will be ignored. |
| `--qmax`           | Set maximum q              | q-points in scattering data larger than this value will be ignored. |

### Settings file
The settings file is a handy tool for setting new default values. All arguments (and even some additional options) can be set in this file, and then automatically be used in future runs. The directory containing the 

The most useful options beyond the arguments from the argument list are:
|  Name  | Description |
|--------|-------------|
| `skip` | Skip the first N entries in the SAXS data file. |
| `format` | Set the format of the outputted figures. |

### q-range
The two arguments `--qmin` and `--qmax` specifies the range of q-values used for the fit. Any data points outside this range in the scattering files will be discarded. 

