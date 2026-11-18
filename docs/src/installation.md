---
title: Installing AUSAXS
description: Install the AUSAXS small-angle scattering calculator on Linux, macOS or Windows with pip install pyausaxs.
journey: install
---

# Installing AUSAXS

Install AUSAXS on any platform by opening a terminal and typing:
```bash
python -m pip install "pyausaxs[all]"
```
This installs both the Python package and the `ausaxs` command, together with the optional plotting and graphical interface dependencies. Check that it works by running:
```bash
ausaxs --version
```

If something goes wrong, see the [troubleshooting](troubleshooting.md) page.

??? question "Terminal says `ausaxs` is not found?"

    The folder pip installed `ausaxs` into is not on your `PATH`. This mostly happens on Windows, and sometimes on macOS with the Python that comes with the system. Add the folder to your `PATH` by running the commands for your platform, then open a new terminal and try `ausaxs --version` again:

    === "Windows"

        In PowerShell:
        ```powershell
        $dir = python -c "from importlib.metadata import files; print(next(f.locate().resolve().parent for f in files('pyausaxs') if f.stem == 'ausaxs'))"
        [Environment]::SetEnvironmentVariable("Path", [Environment]::GetEnvironmentVariable("Path", "User") + ";$dir", "User")
        ```

    === "macOS"

        ```bash
        dir=$(python -c "from importlib.metadata import files; print(next(f.locate().resolve().parent for f in files('pyausaxs') if f.stem == 'ausaxs'))")
        echo "export PATH=\"\$PATH:$dir\"" >> ~/.zshrc
        ```

    === "Linux"

        ```bash
        dir=$(python -c "from importlib.metadata import files; print(next(f.locate().resolve().parent for f in files('pyausaxs') if f.stem == 'ausaxs'))")
        echo "export PATH=\"\$PATH:$dir\"" >> ~/.bashrc
        ```

    The rest of this documentation assumes `ausaxs` is on your `PATH`. If you would rather not change it, `python -m pyausaxs` works in place of `ausaxs` in every command, e.g. `python -m pyausaxs --version`.

## Next steps
- **Use a GPU.** AUSAXS runs only on the CPU out of the box. To use your GPU, you will have to [install the GPU driver](gpu.md) first.
- **Learn by example.** The [tutorials](tutorials.md) walk through complete analyses of real data, from fitting protein structures to validating EM maps.
- **Learn what AUSAXS calculates.** See the pages on [SAXS](saxs.md), [SANS](sans.md), [crystals and nanoparticles](crystals.md) and [powder diffraction](powder-diffraction.md).
- **Look up an option.** Every argument of the command-line tools is listed in the [documentation](saxs-fitter.md).

## Other ways to install
If you cannot use Python, you can instead [download a precompiled executable](precompiled-binaries.md) or [compile AUSAXS from source](compiling.md).
