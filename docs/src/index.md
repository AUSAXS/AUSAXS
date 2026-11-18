---
title: "AUSAXS: the fast small-angle scattering calculator"
description: Documentation for AUSAXS, a fast small-angle scattering calculator built on the Debye equation. Calculate SAXS and SANS profiles from PDB, mmCIF, CIF and XYZ structures, and fit atomic models and EM maps to experimental data.
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

[Get started](getting-started/installation.md){ .md-button .md-button--primary }
[Why AUSAXS?](why-ausaxs.md){ .md-button }

</div>

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> is a fast Debye-equation calculator for small-angle X-ray and neutron scattering. It calculates the [SAXS](applications/saxs.md) and [SANS](applications/sans.md) profiles of proteins and other biomolecules directly from their atomic structure, with an explicit hydration shell and excluded volume, and fits them to experimental data. The same engine calculates the scattering and [powder diffraction](applications/powder-diffraction.md) patterns of [crystals and nanoparticles](applications/crystals.md). It runs on Linux, macOS and Windows, from the command line, from Python, or through a graphical interface.

## Getting started
1. [Install AUSAXS](getting-started/installation.md) with a single command:
    ```bash
    python -m pip install "pyausaxs[all]"
    ```
2. Optionally, [install the GPU driver](getting-started/gpu.md) to run the calculations on a GPU.
3. Follow one of the [tutorials](tutorials.md), which walk through complete analyses of real data.

If something goes wrong, see the [troubleshooting](getting-started/troubleshooting.md) page.

## Tools

<div class="grid cards" markdown>

-   :ausaxs-xray:{ .lg .middle } **[SAXS fitting](documentation/cli/fit.md)**

    ---

    Calculate the SAXS profile of an atomic model and fit it to experimental data.

-   :material-atom:{ .lg .middle } **[SANS fitting](documentation/cli/sans.md)**

    ---

    Calculate the SANS profile of an atomic model and fit it to experimental data, with control over the D₂O fraction of the solvent.

-   :ausaxs-debye:{ .lg .middle } **[Debye scattering](documentation/cli/debye.md)**

    ---

    Calculate the scattering of crystals, nanoparticles and other structures without solvent.

-   :ausaxs-em:{ .lg .middle } **[EM map validation](documentation/cli/em-validation.md)**

    ---

    Validate cryo-EM maps against solution SAXS data using dummy structures built from the map itself.

-   :ausaxs-rigidbody:{ .lg .middle } **[Rigid-body refinement](documentation/cli/rigidbody.md)**

    ---

    Refine multi-domain proteins and complexes against experimental data.

</div>

## Python and utilities

<div class="grid cards" markdown>

-   :material-language-python:{ .lg .middle } **[Python API](https://github.com/AUSAXS/pyAUSAXS)**

    ---

    Call AUSAXS from your own Python scripts and notebooks.

-   :material-monitor-dashboard:{ .lg .middle } **[Graphical interface](documentation/cli/fit.md#graphical-user-interface)**

    ---

    Fit models to SAXS data interactively, without the command line.

-   :material-expansion-card:{ .lg .middle } **[GPU driver](getting-started/gpu.md)**

    ---

    Install and enable the optional GPU driver, and choose which GPU to run on.

-   :material-chart-line:{ .lg .middle } **[CLI plotting](documentation/cli/plotting-scripts.md)**

    ---

    Plot the results of the command-line tools, from fits to $\chi^2$ landscapes.

-   :material-cog-outline:{ .lg .middle } **[Developer setup](getting-started/compiling.md#using-your-build-from-python)**

    ---

    Point the Python package at your own build of the AUSAXS library.

</div>
