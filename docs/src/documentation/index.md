---
title: Documentation - AUSAXS
description: Reference documentation for the AUSAXS command-line tools, the GPU driver, building AUSAXS for development, and how to cite it.
---

# Documentation

This section is the reference for the AUSAXS command-line tools: every argument, output file and configuration option. For an introduction to what each tool is used for, see the [applications](../applications/index.md); for worked examples on real data, see the [tutorials](../tutorials.md).

## Command-line tools

All tools are subcommands of `ausaxs`. Run `ausaxs <command> --help` for a short summary of the arguments.

| Command | Description |
|---------|-------------|
| [`ausaxs fit`](cli/fit.md) | Calculate the SAXS profile of an atomic model and fit it to experimental data. `ausaxs saxs` is an alias. |
| [`ausaxs sans`](cli/sans.md) | Calculate the SANS profile of an atomic model and fit it to experimental data. |
| [`ausaxs debye`](cli/debye.md) | Calculate the Debye scattering of a structure without solvent, such as a crystal or nanoparticle. |
| [`ausaxs em`](cli/em-validation.md) | Validate a cryo-EM map against SAXS data. |
| [`ausaxs rigidbody`](cli/rigidbody.md) | Refine a multi-domain model against SAXS data, driven by a configuration script. |
| [`ausaxs plot`](cli/plotting-scripts.md) | Make figures from the output of the other tools. |
| [`ausaxs gpu`](cli/gpu-driver.md) | Check, enable, disable or remove the GPU driver. |
| [`ausaxs gui`](cli/fit.md#graphical-user-interface) | Start the graphical interface. |

## Developers

- [Compiling from source](../getting-started/compiling.md): build the library yourself, and point the Python package at your build with `ausaxs setup`.
- [macOS debug build](developers/macos-debug-build.md): compile AUSAXS in debug mode on macOS.

## Citing AUSAXS

If you use AUSAXS in published work, see [Citing AUSAXS](citing.md) for the references and BibTeX entries.

## Getting help

If something does not work, check the [troubleshooting](../getting-started/troubleshooting.md) page first. Bugs can be reported on the [issue tracker](https://github.com/AUSAXS/AUSAXS/issues), and questions asked in the [discussions](https://github.com/AUSAXS/AUSAXS/discussions).
