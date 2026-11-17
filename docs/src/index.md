---
title: "AUSAXS: Debye-equation SAXS, WAXS and SANS calculator"
description: Documentation for AUSAXS, a fast small-angle scattering calculator built on the Debye equation. Calculate SAXS, WAXS and SANS profiles from PDB, mmCIF, CIF and XYZ structures, and fit atomic models and EM maps to experimental data.
hide:
  - navigation
  - toc
---

<div class="ausaxs-hero" markdown>

<div class="ausaxs-lockup">
  <!-- Inline copy of assets/logo.svg, so extra.css can theme and animate it. -->
  <svg class="ausaxs-mark" xmlns="http://www.w3.org/2000/svg" viewBox="0 6 60 52" aria-hidden="true">
    <defs>
      <radialGradient id="ausaxs-sphere" cx="0.36" cy="0.32" r="0.75">
        <stop offset="0" stop-color="#fff1b8"/>
        <stop offset="0.45" stop-color="#e3b341"/>
        <stop offset="1" stop-color="#8f620b"/>
      </radialGradient>
    </defs>
    <line class="beam" x1="2" y1="32" x2="10" y2="32" stroke-width="2.6" stroke-linecap="round" opacity="0.45"/>
    <circle class="particle" cx="19" cy="32" r="9.5" fill="url(#ausaxs-sphere)"/>
    <path class="wave wave-1" d="M33.55 19.79 A19 19 0 0 1 33.55 44.21" fill="none" stroke-width="4.2" stroke-linecap="round" opacity="1"/>
    <path class="wave wave-2" d="M40.45 14.00 A28 28 0 0 1 40.45 50.00" fill="none" stroke-width="3.6" stroke-linecap="round" opacity="0.68"/>
    <path class="wave wave-3" d="M47.34 8.22 A37 37 0 0 1 47.34 55.78" fill="none" stroke-width="3.0" stroke-linecap="round" opacity="0.4"/>
  </svg>
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