---
title: Rigid-body refinement against SAXS data - AUSAXS
description: How AUSAXS refines multi-domain proteins and symmetric complexes against SAXS data, by moving rigid bodies of an atomic model under constraints until its profile fits the data.
journey: apps
hide:
  - toc
---

# Rigid-body refinement

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> refines multi-domain proteins and complexes against SAXS or SANS data, by moving rigid parts of an atomic model until its profile fits the measurement.

## Why rigid bodies
- **What is known.** High-resolution structures are often available for the individual domains of a protein, or the subunits of a complex. What is missing is how they are arranged in solution, for instance when they are joined by flexible linkers or crystallised separately.
- **What the data can tell.** A solution scattering profile is a one-dimensional, orientationally averaged measurement. It cannot determine an atomic structure, but it is sensitive to the overall shape. This is what the arrangement of a few large rigid bodies changes.
- **Constraints.** The arrangement must also be physically sensible: domains of the same chain stay connected through their linkers, and bodies may not overlap. Such constraints rule out many arrangements that would otherwise fit equally well.
- **Symmetries.** Many complexes are symmetric assemblies of identical subunits. Building them from a single subunit and refining only the symmetry leaves far fewer parameters to determine from the data.

<figure markdown="span">
![The five rigid bodies of the LAR1-2 protein, drawn as backbone traces in different colours, shift and rotate relative to each other while staying connected as the refinement fits them to SAXS data.](../assets/figures/rigidbody-trajectory.gif#only-light){ width="400" }
![The five rigid bodies of the LAR1-2 protein, drawn as backbone traces in different colours, shift and rotate relative to each other while staying connected as the refinement fits them to SAXS data.](../assets/figures/rigidbody-trajectory-dark.gif#only-dark){ width="400" }
</figure>

Since the data hold only low-resolution information, a good fit does not prove that a model is correct. It should be checked against what else is known about the system.

> [!IMPORTANT]
> Rigid-body refinement is still in active development.

## Example
Here we build the dimer of the SARS-CoV-2 main protease from a single chain, [SASDJG5_single.pdb](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/SASDJG5_single.pdb), by refining a two-fold symmetry against a SAXS measurement of the dimer, [SASDJG5.dat](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/SASDJG5.dat). Save the following script as `dimer.conf`. It declares the symmetry, keeps the two copies within 30 Å of each other, and runs 1000 steps:

> [!TIP]
> We highly recommend the [graphical interface](https://sastutorials.org/RigidBody/RigidBodyGUI.html) for rigid-body refinement. It lets you define the bodies, constraints and symmetries visually and writes the configuration script for you, which makes setting up a refinement much simpler than writing the script by hand.

```
output output/dimer
load {
	pdb SASDJG5_single.pdb
	saxs SASDJG5.dat
}
symmetry p2
constrain {
	first b1s1
	second b1
	type attract
	distance 30
}
save initial.fit

parameters {
	iterations 1000
}
loop
	optimize_once
		on_improvement
			print "{iteration}/{iterations_total}: Accepted with new chi2 {chi2}"
		end
	end
end
save final.fit
save dimer.pdb
```

=== "Command line"

    ```bash
    ausaxs rigidbody dimer.conf
    ```

=== "pyAUSAXS"

    ```python
    import pyausaxs as ausaxs
    with open("dimer.conf") as f:
        rb = ausaxs.prepare_rigidbody_refinement(f.read())
    fit = rb.run() # columns q, I, I_err, I_fit
    ```

Before refinement, the two copies sit too far apart to fit the data; the refined dimer fits it well. Each run follows its own random path, so the final $\chi^2$ varies a little from run to run.

![Before refinement, the two copies of the monomer sit too far apart and fit the data poorly; the refined dimer follows the measured data across the full q range.](../assets/figures/rigidbody-fit.svg#only-light){ width="100%" }
![Before refinement, the two copies of the monomer sit too far apart and fit the data poorly; the refined dimer follows the measured data across the full q range.](../assets/figures/rigidbody-fit-dark.svg#only-dark){ width="100%" }

For a walkthrough in the graphical interface, see the [rigid-body tutorial](https://sastutorials.org/RigidBody/RigidBodyGUI.html). Every option of `ausaxs rigidbody`, and every element of the configuration scripts, is listed in its [documentation](../documentation/cli/rigidbody.md).
