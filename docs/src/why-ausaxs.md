---
title: Why AUSAXS? - AUSAXS
description: Why to use AUSAXS for SAXS and SANS profile calculation. It is fast on CPU and GPU, uses the accurate Debye equation, and reports everything it does, assumes and fits.
hide:
  - toc
---

# Why <span class="ausaxs-name"><span class="au">AU</span>SAXS</span>?

AUSAXS is built on three principles: use the full Debye equation, make it fast, and tell the user everything it does.

## Speed

A complete SAXS fit of a typical protein, hydration shell included, takes a few tens of milliseconds on a desktop CPU. In about 10 seconds, a good consumer CPU handles 300 000 atoms, and a high-end GPU with the optional [GPU driver](getting-started/gpu.md) a million. This makes it practical to fit thousands of models, such as the frames of a molecular dynamics trajectory, or to run long [rigid-body refinements](applications/rigid-body-refinement.md).

=== "SAXS fitting"

    ![A full SAXS fit with AUSAXS takes 16 milliseconds for 600 atoms and 0.6 seconds for 44 000 atoms. Pepsi-SAXS takes about the same time as AUSAXS on small structures and up to four times longer on larger ones; FoXS and CRYSOL are several to over a hundred times slower.](assets/figures/benchmark-scaling-sas.svg#only-light){ width="100%" }
    ![A full SAXS fit with AUSAXS takes 16 milliseconds for 600 atoms and 0.6 seconds for 44 000 atoms. Pepsi-SAXS takes about the same time as AUSAXS on small structures and up to four times longer on larger ones; FoXS and CRYSOL are several to over a hundred times slower.](assets/figures/benchmark-scaling-sas-dark.svg#only-dark){ width="100%" }

=== "Crystals and powder"

    ![The Debye sum of a silver nanoparticle with 800 000 atoms takes AUSAXS 100 seconds on the CPU and 46 seconds on the GPU, against 254 and 48 seconds for AES-Debye. DebyeCalculator is up to 30 times slower than AUSAXS and exceeds the time limit above 200 000 atoms.](assets/figures/benchmark-scaling-debye.svg#only-light){ width="100%" }
    ![The Debye sum of a silver nanoparticle with 800 000 atoms takes AUSAXS 100 seconds on the CPU and 46 seconds on the GPU, against 254 and 48 seconds for AES-Debye. DebyeCalculator is up to 30 times slower than AUSAXS and exceeds the time limit above 200 000 atoms.](assets/figures/benchmark-scaling-debye-dark.svg#only-dark){ width="100%" }

=== "Modern hardware"

    ![AUSAXS alone on a current desktop, with the 2017 desktop for comparison. The CPU calculates the Debye sum of 290 000 atoms in 6 seconds, twice as fast as the 2017 desktop. The GPU, at a bin width of 0.25 Å, takes 1.9 seconds for 830 000 atoms, 15 times faster than the 2017 GPU, and 4.9 seconds for 1.3 million.](assets/figures/benchmark-scaling-modern.svg#only-light){ width="100%" }
    ![AUSAXS alone on a current desktop, with the 2017 desktop for comparison. The CPU calculates the Debye sum of 290 000 atoms in 6 seconds, twice as fast as the 2017 desktop. The GPU, at a bin width of 0.25 Å, takes 1.9 seconds for 830 000 atoms, 15 times faster than the 2017 GPU, and 4.9 seconds for 1.3 million.](assets/figures/benchmark-scaling-modern-dark.svg#only-dark){ width="100%" }

??? info "How the timings were measured"

    Unless noted, timed on a 2017 desktop: an Intel Core i7-8700 (6 cores) and an NVIDIA GeForce GTX 1070 Ti. Each point is the mean of five runs.

    - **Fitting SAXS data:** a full fit of each structure to its SASBDB data, including reading the files and generating the hydration shell, run from the command line with each program's default settings.
    - **Debye sum:** the scattering of spherical silver nanoparticles in vacuum, from the atomic coordinates in memory to I(q) at 1000 points up to q = 10 Å⁻¹. AUSAXS used a bin width of 0.05 Å. An open marker without a line is a single run over the 300-second time limit.
    - **Modern hardware:** AUSAXS alone on a current desktop, an AMD Ryzen 5 9600X (6 cores) and an AMD Radeon RX 9070 XT, with the 2017 desktop in grey. The GPU is shown at a bin width of 0.25 Å. Smaller bin widths may be considerably slower on the GPU for large particles; this is being worked on.
    - **Versions:** Pepsi-SAXS 2020-07-28, FoXS from IMP (develop branch, 11a86f9), CRYSOL from ATSAS 3.2.1, DebyeCalculator 1.1.1 and AES-Debye 1.0.0. The benchmark itself is in [`tests/benchmark/`](https://github.com/AUSAXS/AUSAXS/tree/master/tests/benchmark).

The exact gradient of the Debye sum of point scatterers with respect to the atom positions is also available (`Molecule.debye_raw_vjp` in Python), for example as a loss in machine learning:

![The gradient of the Debye sum with AUSAXS on the CPU takes 5 milliseconds for 225 atoms and 77 seconds for 425 000 atoms. DebyeCalculator takes 37 to 44 times longer at 100 000 atoms on its GPU and CPU and exceeds the time limit at 200 000 atoms. For large structures, the AUSAXS gradient costs about three times as much as calculating I(q).](assets/figures/benchmark-gradient.svg#only-light){ width="100%" }
![The gradient of the Debye sum with AUSAXS on the CPU takes 5 milliseconds for 225 atoms and 77 seconds for 425 000 atoms. DebyeCalculator takes 37 to 44 times longer at 100 000 atoms on its GPU and CPU and exceeds the time limit at 200 000 atoms. For large structures, the AUSAXS gradient costs about three times as much as calculating I(q).](assets/figures/benchmark-gradient-dark.svg#only-dark){ width="100%" }

## Accuracy

Thanks to its weighted bins ([Lytje & Pedersen, 2025](https://doi.org/10.1107/S160057672500562X)), AUSAXS is accurate over a wide range of angles. For SAS, the reference is the exact vacuum profile of the structure, without hydration shell or excluded volume, as these are modelled differently by every program and have no exact solution.

=== "Small angles"

    ![For a protein with 5380 atoms, CRYSOL deviates by 10⁻⁶ from the exact Debye sum up to q = 0.25 inverse ångström and by more than 10 % above q = 0.4. Pepsi-SAXS deviates by 10⁻⁵ to 10⁻³ up to q = 0.65, then by up to 10 % at q = 1. AUSAXS and FoXS deviate by less than 10⁻⁴ up to q = 0.2 and by about 10⁻³ beyond, FoXS rising to a few per cent above q = 0.75.](assets/figures/benchmark-accuracy-sas.svg#only-light){ width="100%" }
    ![For a protein with 5380 atoms, CRYSOL deviates by 10⁻⁶ from the exact Debye sum up to q = 0.25 inverse ångström and by more than 10 % above q = 0.4. Pepsi-SAXS deviates by 10⁻⁵ to 10⁻³ up to q = 0.65, then by up to 10 % at q = 1. AUSAXS and FoXS deviate by less than 10⁻⁴ up to q = 0.2 and by about 10⁻³ beyond, FoXS rising to a few per cent above q = 0.75.](assets/figures/benchmark-accuracy-sas-dark.svg#only-dark){ width="100%" }

=== "Wide angles"

    ![For a silver nanoparticle with 6603 atoms, the relative deviation of AUSAXS from the exact Debye sum lies mostly between 10⁻⁶ and 10⁻⁴ up to q = 10 inverse ångström, almost exactly on top of DebyeCalculator. AES-Debye lies near 10⁻⁸.](assets/figures/benchmark-accuracy-powder.svg#only-light){ width="100%" }
    ![For a silver nanoparticle with 6603 atoms, the relative deviation of AUSAXS from the exact Debye sum lies mostly between 10⁻⁶ and 10⁻⁴ up to q = 10 inverse ångström, almost exactly on top of DebyeCalculator. AES-Debye lies near 10⁻⁸.](assets/figures/benchmark-accuracy-powder-dark.svg#only-dark){ width="100%" }

Beyond a certain point, more accuracy gains nothing: rounding the coordinates to the 0.001 Å of a PDB file already changes I(q) by 10⁻³, and measured SAS data are far less certain still.

??? info "How the accuracy was measured"

    - **Reference:** the exact Debye sum in double precision, computed separately from every program.
    - **Small angles:** the vacuum profile of the protein of SASBDB entry SASDME4 from each SAS program, at its default settings, except that the q range is extended to the program's maximum. As every atom is given the same element, dividing by the profile of a single atom from the same program removes its form factor. AUSAXS uses its default bin width for SAS, 0.25 Å.
    - **Wide angles:** the silver nanoparticle of the speed benchmark. AUSAXS uses a bin width of 0.05 Å.
    - The curves are running medians over 0.05 Å⁻¹ (small angles) and 0.25 Å⁻¹ (wide angles), as the deviation spikes at the minima of I(q).

## Transparency

A fitted profile is only as trustworthy as its assumptions, so AUSAXS tells you what it did, assumed and fitted:

- **Every step and assumption is reported,** such as added implicit hydrogens or a guessed q unit.
- **Unknown input stops the run** instead of being silently guessed; you can override this explicitly.
- **Three fit parameters by default** (scaling, background, hydration shell), all reported with uncertainties, $\chi^2$ and convergence.
- **The model can be inspected:** the structure with its hydration shell is written out, and the profile is split into atom and hydration contributions.
- **Settings can be saved,** edited and reused.

And it is [open source](https://github.com/AUSAXS/AUSAXS), so anything else can be checked in the code.

## One engine, many applications

The same calculation drives [SAXS](applications/saxs.md), [SANS](applications/sans.md), [rigid-body refinement](applications/rigid-body-refinement.md), [EM map validation](applications/em.md), [crystals and nanoparticles](applications/crystals.md) and [powder diffraction](applications/powder-diffraction.md). It runs on Linux, macOS and Windows, from the command line, [Python](https://github.com/AUSAXS/pyAUSAXS) or a [graphical interface](documentation/cli/fit.md#graphical-user-interface), and [installs](getting-started/installation.md) with a single command.
