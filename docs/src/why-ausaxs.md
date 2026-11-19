---
title: Why AUSAXS? - AUSAXS
description: Why to use AUSAXS for SAXS and SANS profile calculation. It is fast on CPU and GPU, uses the full Debye equation, and reports everything it does, assumes and fits.
hide:
  - toc
---

# Why <span class="ausaxs-name"><span class="au">AU</span>SAXS</span>?

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> is built on three principles: use the full Debye equation, make it fast, and tell the user everything we do.

## Speed

A complete SAXS fit, including generating the hydration shell, typically takes a few tens of milliseconds for a protein on an ordinary desktop CPU. Every calculation is multithreaded, so even very large structures remain practical: a good consumer CPU handles about 300 000 atoms in 10 seconds. The optional [GPU driver](getting-started/gpu.md) goes further still: a high-end consumer GPU handles a million atoms in the same time. This makes it practical to fit thousands of models, such as frames of a molecular dynamics trajectory, or to run a [rigid-body refinement](applications/rigid-body-refinement.md) with many thousands of steps.

=== "All"

    ![On a log-log plot of time against the number of atoms, AUSAXS fits SAXS data about as fast as Pepsi-SAXS on small structures and faster on large ones, and several to over a hundred times faster than FoXS and CRYSOL. For the Debye sum of nanoparticles, from 200 to 800 000 atoms, AUSAXS is faster than DebyeCalculator and AES-Debye on the CPU at every size, and on the GPU it matches AES-Debye at the largest sizes.](assets/figures/benchmark-scaling.svg#only-light){ width="100%" }
    ![On a log-log plot of time against the number of atoms, AUSAXS fits SAXS data about as fast as Pepsi-SAXS on small structures and faster on large ones, and several to over a hundred times faster than FoXS and CRYSOL. For the Debye sum of nanoparticles, from 200 to 800 000 atoms, AUSAXS is faster than DebyeCalculator and AES-Debye on the CPU at every size, and on the GPU it matches AES-Debye at the largest sizes.](assets/figures/benchmark-scaling-dark.svg#only-dark){ width="100%" }

=== "SAXS fitting"

    ![A full SAXS fit with AUSAXS takes 16 milliseconds for 600 atoms and 0.6 seconds for 44 000 atoms. Pepsi-SAXS takes about the same time as AUSAXS on small structures and up to four times longer on larger ones; FoXS and CRYSOL are several to over a hundred times slower.](assets/figures/benchmark-scaling-sas.svg#only-light){ width="100%" }
    ![A full SAXS fit with AUSAXS takes 16 milliseconds for 600 atoms and 0.6 seconds for 44 000 atoms. Pepsi-SAXS takes about the same time as AUSAXS on small structures and up to four times longer on larger ones; FoXS and CRYSOL are several to over a hundred times slower.](assets/figures/benchmark-scaling-sas-dark.svg#only-dark){ width="100%" }

=== "Crystals and powder"

    ![The Debye sum of a silver nanoparticle with 800 000 atoms takes AUSAXS 100 seconds on the CPU and 46 seconds on the GPU, against 254 and 48 seconds for AES-Debye. DebyeCalculator is up to 30 times slower than AUSAXS and exceeds the time limit above 200 000 atoms.](assets/figures/benchmark-scaling-debye.svg#only-light){ width="100%" }
    ![The Debye sum of a silver nanoparticle with 800 000 atoms takes AUSAXS 100 seconds on the CPU and 46 seconds on the GPU, against 254 and 48 seconds for AES-Debye. DebyeCalculator is up to 30 times slower than AUSAXS and exceeds the time limit above 200 000 atoms.](assets/figures/benchmark-scaling-debye-dark.svg#only-dark){ width="100%" }

=== "Modern hardware"

    ![AUSAXS alone on a 2024 desktop, with the 2017 desktop for comparison. The CPU calculates the Debye sum of 290 000 atoms in 6 seconds, twice as fast as the 2017 desktop. The GPU, at a bin width of 0.25 Å, takes 1.9 seconds for 830 000 atoms, 15 times faster than the 2017 GPU, and 4.9 seconds for 1.3 million.](assets/figures/benchmark-scaling-modern.svg#only-light){ width="100%" }
    ![AUSAXS alone on a 2024 desktop, with the 2017 desktop for comparison. The CPU calculates the Debye sum of 290 000 atoms in 6 seconds, twice as fast as the 2017 desktop. The GPU, at a bin width of 0.25 Å, takes 1.9 seconds for 830 000 atoms, 15 times faster than the 2017 GPU, and 4.9 seconds for 1.3 million.](assets/figures/benchmark-scaling-modern-dark.svg#only-dark){ width="100%" }

Unless noted, timed on a 2017 desktop: an Intel Core i7-8700 (6 cores) and an NVIDIA GeForce GTX 1070 Ti. Each point is the mean of five runs.

- **Fitting SAXS data:** a full fit of each structure to its SASBDB data, including reading the files and generating the hydration shell, run from the command line with each program's default settings. By default, AUSAXS fits q up to 0.5 Å⁻¹, while the other programs fit the whole data file.
- **Debye sum:** the scattering of spherical silver nanoparticles in vacuum, from the atomic coordinates in memory to I(q) at 1000 points up to q = 10 Å⁻¹. AUSAXS used its default bin width for crystals, 0.05 Å. An open marker without a line is a single run over the 300-second time limit.
- **Modern hardware:** <span class="ausaxs-name"><span class="au">AU</span>SAXS</span> alone on a 2024 desktop, an AMD Ryzen 5 9600X (6 cores) and an AMD Radeon RX 9070 XT, with the 2017 desktop in grey. The GPU is shown at a bin width of 0.25 Å. Smaller bin widths may be considerably slower on the GPU for large particles; this is being worked on.
- **Versions:** Pepsi-SAXS 2020-07-28, FoXS from IMP (develop branch, 11a86f9), CRYSOL from ATSAS 3.2.1, DebyeCalculator 1.1.1 and AES-Debye 1.0.0. The benchmark itself is in [`tests/benchmark/`](https://github.com/AUSAXS/AUSAXS/tree/master/tests/benchmark).

The scattering curve can also be differentiated with respect to the atom positions, for example to use it as a loss in machine learning. In Python this is `Molecule.debye_raw_vjp`. The gradient is that of the exact Debye sum, and on the CPU it is up to 40 times faster than DebyeCalculator on either CPU or GPU:

![The gradient of the Debye sum with AUSAXS on the CPU takes 5 milliseconds for 225 atoms and 77 seconds for 425 000 atoms. DebyeCalculator takes 37 to 44 times longer at 100 000 atoms on its GPU and CPU and exceeds the time limit at 200 000 atoms. For large structures, the AUSAXS gradient costs about three times as much as calculating I(q).](assets/figures/benchmark-gradient.svg#only-light){ width="100%" }
![The gradient of the Debye sum with AUSAXS on the CPU takes 5 milliseconds for 225 atoms and 77 seconds for 425 000 atoms. DebyeCalculator takes 37 to 44 times longer at 100 000 atoms on its GPU and CPU and exceeds the time limit at 200 000 atoms. For large structures, the AUSAXS gradient costs about three times as much as calculating I(q).](assets/figures/benchmark-gradient-dark.svg#only-dark){ width="100%" }

## The full Debye equation

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> sums the contribution of every pair of atoms, rather than relying on a multipole expansion or on coarse-graining the structure. The hydration shell and the excluded volume are modelled explicitly in real space. How this is done, and how it compares with other programs, is described in [Lytje & Pedersen (2025)](https://doi.org/10.1107/S160057672500562X).

To be fast, the pair distances are sorted into bins before the sum over q. The error this introduces shrinks with the square of the bin width. At the default settings for SAS, it is below 10⁻⁴ up to q = 0.25 Å⁻¹ and about 10⁻³ beyond. Programs based on a multipole expansion, such as CRYSOL and Pepsi-SAXS, are very accurate at small angles, but fail abruptly once the expansion order no longer suffices. At the finer default bins for crystals and powder diffraction, the deviation is about 10⁻⁵ over the whole range up to q = 10 Å⁻¹, the same as DebyeCalculator:

=== "All"

    ![The relative deviation from the exact Debye sum of every program, against q from 0.01 to 10 inverse ångström on a logarithmic axis. For the protein, CRYSOL is accurate to 10⁻⁶ but fails above q = 0.3, Pepsi-SAXS lies between 10⁻⁵ and 10⁻³ but fails above q = 0.65, and FoXS and AUSAXS rise to about 10⁻³ at q = 0.3 and stay within 10⁻² up to q = 2. For the silver nanoparticle, AUSAXS and DebyeCalculator lie between 10⁻⁶ and 10⁻⁴ up to q = 10, and AES-Debye near 10⁻⁸.](assets/figures/benchmark-accuracy.svg#only-light){ width="100%" }
    ![The relative deviation from the exact Debye sum of every program, against q from 0.01 to 10 inverse ångström on a logarithmic axis. For the protein, CRYSOL is accurate to 10⁻⁶ but fails above q = 0.3, Pepsi-SAXS lies between 10⁻⁵ and 10⁻³ but fails above q = 0.65, and FoXS and AUSAXS rise to about 10⁻³ at q = 0.3 and stay within 10⁻² up to q = 2. For the silver nanoparticle, AUSAXS and DebyeCalculator lie between 10⁻⁶ and 10⁻⁴ up to q = 10, and AES-Debye near 10⁻⁸.](assets/figures/benchmark-accuracy-dark.svg#only-dark){ width="100%" }

=== "Small angles"

    ![For a protein with 5380 atoms, CRYSOL deviates by 10⁻⁶ from the exact Debye sum up to q = 0.25 inverse ångström and by more than 10 % above q = 0.4. Pepsi-SAXS deviates by 10⁻⁵ to 10⁻³ up to q = 0.65, then by up to 10 % at q = 1. AUSAXS and FoXS deviate by less than 10⁻⁴ up to q = 0.2 and by about 10⁻³ beyond, FoXS rising to a few per cent above q = 0.75.](assets/figures/benchmark-accuracy-sas.svg#only-light){ width="100%" }
    ![For a protein with 5380 atoms, CRYSOL deviates by 10⁻⁶ from the exact Debye sum up to q = 0.25 inverse ångström and by more than 10 % above q = 0.4. Pepsi-SAXS deviates by 10⁻⁵ to 10⁻³ up to q = 0.65, then by up to 10 % at q = 1. AUSAXS and FoXS deviate by less than 10⁻⁴ up to q = 0.2 and by about 10⁻³ beyond, FoXS rising to a few per cent above q = 0.75.](assets/figures/benchmark-accuracy-sas-dark.svg#only-dark){ width="100%" }

=== "Wide angles"

    ![For a silver nanoparticle with 6603 atoms, the relative deviation of AUSAXS from the exact Debye sum lies mostly between 10⁻⁶ and 10⁻⁴ up to q = 10 inverse ångström, almost exactly on top of DebyeCalculator. AES-Debye lies near 10⁻⁸.](assets/figures/benchmark-accuracy-powder.svg#only-light){ width="100%" }
    ![For a silver nanoparticle with 6603 atoms, the relative deviation of AUSAXS from the exact Debye sum lies mostly between 10⁻⁶ and 10⁻⁴ up to q = 10 inverse ångström, almost exactly on top of DebyeCalculator. AES-Debye lies near 10⁻⁸.](assets/figures/benchmark-accuracy-powder-dark.svg#only-dark){ width="100%" }

- **Reference:** the exact Debye sum in double precision, computed separately from every program.
- **Small angles:** the vacuum profile of the protein SASDME4 from each SAS program, at its default settings, except that the q range is extended to the program's maximum. As every atom is given the same element, the profile of a single atom from the same program removes its form factor. AUSAXS uses its default bin width for SAS, 0.25 Å.
- **Wide angles:** the silver nanoparticle of the speed benchmark. AUSAXS uses its default bin width for crystals, 0.05 Å.
- The curves are running medians over 0.05 Å⁻¹ (small angles) and 0.25 Å⁻¹ (wide angles), as the deviation spikes at the minima of I(q).

Beyond a certain point, a more accurate sum gains nothing. Merely writing a structure to a PDB file, which rounds the coordinates to 0.001 Å, already changes I(q) by 10⁻³ or more, so a precision much better than 10⁻⁴ has no physical meaning. For SAS, even an error of 1–2 % is acceptable, as measured data are always considerably less certain. Moreover, the Debye sum of independent atoms is itself only a model. It leaves out effects that matter more and more at wide angles, so the true error of any calculated profile there is larger than these curves suggest.

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
