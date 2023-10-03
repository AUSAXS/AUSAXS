# Validation of electron microscopy maps
Welcome to this short tutorial on our EM validation method. Here we will cover all of the arguments and settings which are available for the fitting procedure. 

The program requires just two mandatory arguments: a location of a map file, and a location of a SAXS data file. 
`em_fitter <map file> <saxs file>`

## Output
By default all output will be saved in a new `output/em_fitter/<map-name>/` folder, although this location can be changed by specifying the `--output` argument (note that the subfolder `<map-name>` will always be created). 

In this folder you will find a bunch of `.plot` files. These are instructions for making plots with the `plot.py` python script. Simply run the plotting script with this folder as its argument. You can optionally add the `--big` (experimental) argument for larger fonts. After you have plotted the files, you will have the three following $\chi^2$-landscape figures:
- **chi2_evaluated_points_full**: This figure shows the full $\chi^2$-landscape for the full range of scanned threshold values. 
- **chi2_evaluated_points_limited**: This figure shows the same range as the first, but with very high-$\chi^2$ points removed. This is done so the more interesting low-$\chi^2$ range can better be visualized. A rolling average is performed, as shown by the solid red line. This average is also the basis for checking for the presence of more local minima. All minima are shown with a dashed vertical red line. You can find dummy models for each of these minima in the `models` folder.
- **chi2_near_minimum**: This figure shows the area near the absolute minimum, sampled with the smallest step-size possible for your map. Shown are also the mean and standard deviation of the entire range. This mean is the reported average $\chi^2$ value, $\bar{\chi}^2$. The blue dot is just the interpolated $\chi^2$ value from earlier, and may thus often be somewhat higher/lower in value than the rest of the plot. 

The first and third of these plots also have an accompanying `_mass` version, where the threshold x-axis has been replaced with a mass-axis. 

Beyond these landscape plots, you will also find the four figures:
- **intensity_fit**: A simple plot of the original scattering data with the fit superimposed.
- **residuals**: The fit residuals. 
- **log**: A nicer figure combining the fit and residuals into one, with a linear x-axis. 
- **loglog**: Same as before, but with a logarithmic x-axis. 

You will also find a copy of the *used* parts of the original scattering data file, along with a `fit.fit` file. Both files are compatible with the plotting utility from the `ATSAS` package. The other text file is the `report.txt`, which contains information about the fit itself. Here you will find the actual parameter values, along with the fitted threshold value and the $\chi^2$. 

Finally you will also find the `model.pdb` file, which contains the dummy structure. Note that if your map is very big, this file may be split into several numerated parts. 

## Arguments
Arguments are handled by [CLI11](https://github.com/CLIUtils/CLI11). You can see a list of available arguments by providing the `--help` argument. 
| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`-h`, `--help`      | Help                       | Print a help message |
|`-s`, `--settings`  | Set settings file          | Set a settings file. This file can be used to easily set new default values for most agruments. |
|`-o`, `--output`    | Set output directory       | Set the path to the output directory |
|`-t`, `--threads`   | Set the number of threads  | Use this to limit the number of threads. By default all but one thread available on the system will be used. |
|`--quiet`           | Disbale optional information | No optional progress information will be printed to the console. |
| `--qmin`           | Set minimum q              | q-points in scattering data less than this value will be ignored. |
| `--qmax`           | Set maximum q              | q-points in scattering data larger than this value will be ignored. |
| `--levelmin`       | Set minimum RMS-level      | This determines the minimum threshold level used in the fit. |
| `--levelmax`       | Set maximum RMS-level      | This determines the maximum threshold level used in the fit. |
| `--frequency`      | Set the sampling frequency | This determines how often the map will be sampled for dummy atoms. |
| `--max-iterations` | Set the maximum number of evaluations | This is a soft limit on the number of evaluations made during the fit. Typically up to 50% more evaluations than this value are used. |
| `--no-hydrate`     | Disable hydration          | The dummy models will not be hydrated. |
| `--fixed-weight`   | Enable fixed weights       | The weight of all dummy atoms will be set to unity. |

### Settings file
The settings file is a handy tool for setting new default values. All arguments (and even some additional options) can be set in this file, and then automatically be used in future runs. The settings file can either be manually set using the `--settings` argument, or by placing it in the same folder as the structure file, in which case it will be used automatically. 

The most useful options beyond the arguments from the argument list are 
|  Name  | Description |
|--------|-------------|
| `skip` | Skip the first N entries in the SAXS data file. |
| `charge-levels` | Set the number of partial histograms to use. These will be uniformly distributed between the `levelmin` and `levelmax` levels. This is primarily useful if the scanned range is very large. |
| `format` | Set the format of the outputted figures. |

### q-range
The two arguments `--qmin` and `--qmax` specifies the range of q-values used for the fit. Any data points outside this range in the scattering files will be discarded. 

### level-range
The two arguments `--levelmin` and `--levelmax` specifies the range of threshold values in terms of the RMS. This is the same levels as those used by e.g. PyMOL, meaning the user can use PyMOL for estimating the bounds, and then apply them with these two arguments. Note that the evaluation of higher threshold values are significantly faster than smaller values, meaning the choice of `--levelmin` has a way higher impact on the performance than `--levelmax`. 

### Map sampling frequency
The argument `--frequency` determines how often dummy atoms will be sampled from the grid intrinsic in the map. A value of `2` means that only every second grid point along each axis is used, effectively reducing the size of the dummy model by a factor 8. For somewhat larger high-resolution maps, using `--frequency 2` may be enough for a decent quick fit. 

Note that since this sampling frequency directly relates to the information content in the dummy model, using higher values introduces more noise to the $\chi^2$-landscape. 

### To hydrate or not to hydrate
By default all dummy structures are hydrated. While our hydration algorithm runs in linear time, it *does* have a significant impact on the runtime. If you believe your structure does not need a hydration layer for some reason, this feature can be disabled by adding the `--no-hydrate`

### Using fixed weights
By default all dummy atoms are weighted with the charge density at their specific location in the EM map. This works perfectly well in most cases, but it does have a caveat: low-density noise can easily be included without significantly impacting the evaluated scattering curve. Thus we have added this option which sets all dummy atom weights to unity, meaning the low-density noise will now become relevant again. This can sometimes help with shifting minima from the low-threshold noisy region to more physical higher threshold regions.  



