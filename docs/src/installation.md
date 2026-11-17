---
title: Installing AUSAXS
description: Install the AUSAXS small-angle scattering calculator on Linux, macOS or Windows with pip install pyausaxs.
---

# Installing AUSAXS

AUSAXS is installed with a single command, on Linux, macOS and Windows alike. Open a terminal and run:
```bash
python -m pip install "pyausaxs[all]"
```
This installs both the Python package and the `ausaxs` command, together with the optional plotting and graphical interface dependencies. Check that it works by running:
```bash
ausaxs --version
```

> [!NOTE]
> If your shell cannot find `ausaxs`, the folder pip installs scripts into is not on your path. Use `python -m pyausaxs` in its place, e.g. `python -m pyausaxs --version`.

If something goes wrong, see the [troubleshooting](troubleshooting.md) page.

## Next steps
- **Use a GPU.** AUSAXS runs only on the CPU out of the box. To use your GPU, you will have to [install the GPU backend](gpu.md) first.
- **Learn by example.** The [tutorials](tutorials.md) walk through complete analyses of real data, from fitting protein structures to validating EM maps.
- **Learn what AUSAXS calculates.** See the pages on [SAXS](saxs.md), [SANS](sans.md), [crystals and nanoparticles](crystals.md) and [powder diffraction](powder-diffraction.md).
- **Look up an option.** Every argument of the command-line tools is listed in the [documentation](saxs-fitter.md).

## Other ways to install
If you cannot use Python, you can instead [download a precompiled executable](precompiled-binaries.md) or [compile AUSAXS from source](compiling.md).
