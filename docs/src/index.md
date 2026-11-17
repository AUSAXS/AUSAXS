---
title: "AUSAXS: Debye-equation SAXS, WAXS and SANS calculator"
description: Documentation for AUSAXS, a fast small-angle scattering calculator built on the Debye equation. Calculate SAXS, WAXS and SANS profiles from PDB, mmCIF, CIF and XYZ structures, and fit atomic models and EM maps to experimental data.
hide:
  - navigation
  - toc
---

<div class="ausaxs-hero" markdown>

<div class="ausaxs-lockup">
  <img class="ausaxs-mark" src="assets/logo-light.svg#only-light" alt="">
  <img class="ausaxs-mark" src="assets/logo.svg#only-dark" alt="">
  <h1 class="ausaxs-wordmark"><span class="au">AU</span>SAXS</h1>
</div>

<p class="ausaxs-tagline">Debye scattering from any atomic structure, all in one package.</p>

[Get started](installation.md){ .md-button .md-button--primary }
[Tutorials](tutorials.md){ .md-button }

<div class="ausaxs-pills" markdown="span">[SAXS](saxs.md) [SANS](sans.md) [Crystals](crystals.md) [Powder diffraction](powder-diffraction.md) [Rigid-body refinement](rigidbody.md) [GPU](gpu.md) [Linux · macOS · Windows](installation.md)</div>

</div>

AUSAXS accurately predicts small-angle X-ray scattering (SAXS), wide-angle X-ray scattering (WAXS) and small-angle neutron scattering (SANS) profiles of proteins and other biomolecules, and fits them to experimental small-angle scattering data. The same kernel computes the Debye scattering of crystals, nanoparticles and arbitrary point scatterers, with or without atomic form factors. It is available both as a command-line program and as a Python package:
```bash
python -m pip install pyausaxs
```

## Getting started
1. [Install AUSAXS](installation.md) with a single `pip` command.
2. Optionally, [install the GPU backend](gpu.md) to run the calculations on a GPU.
3. Follow one of the [tutorials](tutorials.md), which walk through complete analyses of real data.

If you have issues running the programs, please take a look at the [troubleshooting](troubleshooting.md) page.

## Tools
All tools are run through the single `ausaxs` executable as `ausaxs <tool> [options]`. Each card links to the documentation of every option of the tool.

<div class="grid cards ausaxs-two-col" markdown>

-   :material-chart-bell-curve-cumulative:{ .lg .middle } **SAXS fitting** &nbsp; `ausaxs fit`

    ---

    Calculate the SAXS profile of a protein or other macromolecule from its atomic model, and test whether it is consistent with experimental SAXS data. The hydration shell and excluded volume are modelled explicitly from the atomic coordinates.

    [:octicons-arrow-right-24: Documentation](saxs-fitter.md)

-   :material-cube-scan:{ .lg .middle } **EM map validation** &nbsp; `ausaxs em`

    ---

    Validate cryo-EM maps against solution SAXS data using dummy structures built from the map itself.

    [:octicons-arrow-right-24: Documentation](em-validation.md)

-   :material-puzzle-outline:{ .lg .middle } **Rigid-body refinement** &nbsp; `ausaxs rigidbody`

    ---

    Refine multi-domain proteins and complexes against SAXS data, driven by configuration scripts. This tool is still in active development.

    [:octicons-arrow-right-24: Documentation](rigidbody.md)

-   :material-language-python:{ .lg .middle } **Python API and GUI**

    ---

    All three tools are also available through a graphical interface and a Python API in the [pyAUSAXS](https://github.com/AUSAXS/pyAUSAXS) package. Its [examples](https://github.com/AUSAXS/pyAUSAXS/tree/master/examples) cover everything from a basic Debye calculation to model fitting and crystal scattering.

    [:octicons-arrow-right-24: pyAUSAXS on GitHub](https://github.com/AUSAXS/pyAUSAXS)

</div>

The [plotting scripts](plotting-scripts.md) turn the output of the command-line tools into figures.

## Citing AUSAXS
If you use AUSAXS in your work, please cite:

- K. Lytje and J. S. Pedersen. Small-angle X-ray scattering profile calculation for high-resolution models of biomacromolecules. _J. Appl. Cryst._ (2025). doi: [10.1107/S160057672500562X](https://doi.org/10.1107/S160057672500562X)

For EM map validation, please also cite:

- Validation of electron-microscopy maps using solution small-angle X-ray scattering. _Acta Cryst. D_ (2024). doi: [10.1107/S2059798324005497](https://doi.org/10.1107/S2059798324005497)
