---
title: Powder diffraction patterns of nanocrystals - AUSAXS
description: Calculate X-ray powder diffraction patterns of nanocrystals and nanoparticles from their atomic structure with the Debye scattering equation in AUSAXS.
---

# Powder diffraction

A powder averages over every orientation of its crystallites, which is exactly what the Debye scattering equation calculates. For a crystallite small enough to list every atom, typically a nanocrystal or nanoparticle, AUSAXS calculates its powder pattern directly from the atomic structure, with no assumptions about peak shapes or long-range order. Larger crystallites are better handled by Rietveld refinement.

## What AUSAXS calculates
The intensity $I(q)$ of the particle, with the X-ray form factor of each element, up to $q = 10$ Å⁻¹. This covers the full range of a laboratory Cu Kα instrument, which reaches about $8$ Å⁻¹.

How to set up the calculation, and how accurate it is, is described on the [crystals and nanoparticles](crystals.md) page.

## Running it
```bash
ausaxs debye particle.cif --qmax 8 --wavelength 1.5406
```
This writes the pattern both as a function of $q$ and of the scattering angle $2\theta$, using

$$
q = \frac{4\pi}{\lambda}\sin\theta .
$$

The options of `ausaxs debye` are listed on the [crystals and nanoparticles](crystals.md#running-it) page.

## What you add yourself
AUSAXS calculates the scattering of the particle alone. To compare it with a measured pattern, apply the corrections of your instrument, such as polarisation and peak broadening.
