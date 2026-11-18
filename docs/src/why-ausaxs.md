---
title: Why AUSAXS? - AUSAXS
description: Why to use AUSAXS for SAXS and SANS profile calculation. It is fast on CPU and GPU, uses the full Debye equation, and reports everything it does, assumes and fits.
hide:
  - toc
---

# Why <span class="ausaxs-name"><span class="au">AU</span>SAXS</span>?

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> is built on three principles: use the full Debye equation, make it fast, and tell the user everything we do.

## Speed

<!--
  Benchmark figure goes here once the paper benchmarks have been run:
  ![...](assets/figures/benchmark-sas.svg#only-light){ width="100%" }
  ![...](assets/figures/benchmark-sas-dark.svg#only-dark){ width="100%" }
  Caption: what is timed (a full fit, including the hydration shell), the hardware, the program versions,
  and a link to tests/benchmark/paper/ so the numbers can be reproduced.
-->

A complete SAXS fit, including generating the hydration shell, typically takes a few tens of milliseconds for a protein on an ordinary desktop CPU. Every calculation is multithreaded, so even very large structures remain practical: a good consumer CPU handles about 300 000 atoms in 10 seconds. The optional [GPU driver](getting-started/gpu.md) goes further still: a high-end consumer GPU handles a million atoms in the same time. This makes it practical to fit thousands of models, such as frames of a molecular dynamics trajectory, or to run a [rigid-body refinement](applications/rigid-body-refinement.md) with many thousands of steps.

## The full Debye equation

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> sums the contribution of every pair of atoms, rather than relying on a multipole expansion or on coarse-graining the structure. The hydration shell and the excluded volume are modelled explicitly in real space. How this is done, and how it compares with other programs, is described in [Lytje & Pedersen (2025)](https://doi.org/10.1107/S160057672500562X).

## Transparency

A fitted profile is only as trustworthy as the assumptions behind it. <span class="ausaxs-name"><span class="au">AU</span>SAXS</span> therefore tells you what it did, what it assumed and what it fitted:

- **Every step is reported.** For example, when implicit hydrogens are added, when hydration atoms already in the file are discarded, or how the columns of a data file were read.
- **Assumptions come with a warning.** If the unit of the q values cannot be determined from the data file, the program says which one it assumed. Automatic choices, such as the binning used for a structure, are also printed.
- **Unknown input stops the program.** An atom or residue that the program does not recognise ends the run, rather than being silently assigned a guessed charge. You can override this, but you have to do so explicitly.
- **Few fit parameters, all reported.** By default, a SAXS fit has three free parameters: the scaling, a constant background and the scattering power of the hydration shell. All of them are reported with their uncertainties, together with $\chi^2$, the number of degrees of freedom and whether the fit converged. Additional parameters are only added when you ask for them.
- **The model can be inspected.** The structure is written out as it was used in the fit, including the generated hydration shell, so you can look at it in any molecular viewer. The profile and distance distribution are split into the contributions from the atoms and from the hydration shell, and the calculated volume, charge and density of the molecule are printed for checking.
- **Settings can be saved.** All settings of a run can be written to a file, edited and reused. The excluded volume and the grid can be saved as structures too.

Finally, <span class="ausaxs-name"><span class="au">AU</span>SAXS</span> is [open source](https://github.com/AUSAXS/AUSAXS), so anything this page does not answer can be checked in the code.

## One engine for many applications

The same calculation is used for [SAXS](applications/saxs.md), [SANS](applications/sans.md), [rigid-body refinement](applications/rigid-body-refinement.md) and [EM map validation](applications/em.md), and without the solvent for [crystals and nanoparticles](applications/crystals.md) and [powder diffraction](applications/powder-diffraction.md). Improvements to the engine therefore benefit every application.

## Use it from anywhere

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> runs on Linux, macOS and Windows. It can be used from the command line, from [Python](https://github.com/AUSAXS/pyAUSAXS), or through a [graphical interface](documentation/cli/fit.md#graphical-user-interface), and it is [installed](getting-started/installation.md) with a single command.
