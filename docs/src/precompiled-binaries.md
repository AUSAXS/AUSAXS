---
title: Precompiled executables - AUSAXS
description: Download a precompiled AUSAXS executable for Linux, macOS or Windows from the GitHub releases, without installing Python.
---

# Precompiled executables

If you cannot use Python, precompiled executables for Linux, macOS and Windows are available as [GitHub releases](https://github.com/AUSAXS/AUSAXS/releases). The [Python package](installation.md) remains the recommended way to install AUSAXS.

1. Download the zip file for your platform from the [latest release](https://github.com/AUSAXS/AUSAXS/releases/latest).
2. Extract it. The `ausaxs` executable is in the `bin` folder.
3. Check that it works by running `ausaxs --version` from that folder in a terminal.

Without the Python package, the `.plot` files written by the tools must be converted to images with the `scripts/plot.py` script from the AUSAXS repository, see [CLI plotting](plotting-scripts.md).

If the executable does not start, see the [troubleshooting](troubleshooting.md) page. Then continue with the [next steps](installation.md#next-steps).
