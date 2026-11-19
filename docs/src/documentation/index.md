---
title: Documentation - AUSAXS
description: Reference documentation for the AUSAXS command-line tools, the GPU driver, building AUSAXS for development, and how to cite it.
---

# Documentation

This section is the reference for the AUSAXS command-line tools: every argument, output file and configuration option. For an introduction to what each tool is used for, see the [applications](../applications/index.md); for worked examples on real data, see the [tutorials](../tutorials.md).

## Command-line tools

All tools are subcommands of `ausaxs`. Run `ausaxs <command> --help` for a short summary of the arguments.

<div class="ausaxs-cmds" markdown>
<div class="ausaxs-cmds__head"><span>Command</span><span>Description</span></div>

[<code>ausaxs fit</code><span class="ausaxs-cmd__desc">Fit an atomic model to SAXS data. Also available as <code>ausaxs saxs</code>.</span>](cli/fit.md){ .ausaxs-cmd }

[<code>ausaxs sans</code><span class="ausaxs-cmd__desc">Fit an atomic model to SANS data.</span>](cli/sans.md){ .ausaxs-cmd }

[<code>ausaxs debye</code><span class="ausaxs-cmd__desc">Calculate the scattering of crystals and nanoparticles, without solvent.</span>](cli/debye.md){ .ausaxs-cmd }

[<code>ausaxs em</code><span class="ausaxs-cmd__desc">Validate a cryo-EM map against SAXS data.</span>](cli/em-validation.md){ .ausaxs-cmd }

[<code>ausaxs rigidbody</code><span class="ausaxs-cmd__desc">Refine the arrangement of domains against SAXS data.</span>](cli/rigidbody.md){ .ausaxs-cmd }

[<code>ausaxs plot</code><span class="ausaxs-cmd__desc">Make figures from the output of the other tools.</span>](cli/plotting-scripts.md){ .ausaxs-cmd }

[<code>ausaxs gpu</code><span class="ausaxs-cmd__desc">Install and manage the GPU driver.</span>](cli/gpu-driver.md){ .ausaxs-cmd }

[<code>ausaxs gui</code><span class="ausaxs-cmd__desc">Start the graphical interface.</span>](cli/fit.md#graphical-user-interface){ .ausaxs-cmd }

[<code>ausaxs setup</code><span class="ausaxs-cmd__desc">Point the Python package at a library you compiled yourself.</span>](../getting-started/compiling.md#using-your-build-from-python){ .ausaxs-cmd }

</div>

> [!TIP]
> If your terminal says `ausaxs` is not found, type `python -m pyausaxs` in its place, e.g. `python -m pyausaxs fit --help`. To use the short `ausaxs` instead, add it to your `PATH` as described on the [installation](../getting-started/installation.md) page.

## Developers

- [Compiling from source](../getting-started/compiling.md): build the library yourself, and point the Python package at your build with `ausaxs setup`.
- [Running the tests](developers/tests.md): build and run the test suite.
- [Contributing](https://github.com/AUSAXS/AUSAXS/blob/master/CONTRIBUTING.md): how to report bugs, make suggestions and contribute code.

## Citing AUSAXS

If you use AUSAXS in published work, see [Citing AUSAXS](citing.md) for the references and BibTeX entries.

## Getting help

If something does not work, check the [troubleshooting](../getting-started/troubleshooting.md) page first. Bugs can be reported on the [issue tracker](https://github.com/AUSAXS/AUSAXS/issues), and questions asked in the [discussions](https://github.com/AUSAXS/AUSAXS/discussions).
