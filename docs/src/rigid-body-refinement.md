---
title: Rigid-body refinement against SAXS data - AUSAXS
description: How AUSAXS refines multi-domain proteins and symmetric complexes against SAXS data, by moving rigid bodies of an atomic model under constraints until its profile fits the data.
journey: apps
hide:
  - toc
---

# Rigid-body refinement

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> refines multi-domain proteins and complexes against SAXS data. Rigid parts of an atomic model are moved around until its profile fits the measurement.

## How it works
- **Bodies.** The structure is split into rigid bodies: by chain, at given residues, or one per file.
- **Constraints.** The bodies are held together by constraints, such as the backbone bonds between the domains of a single chain. Overlapping bodies are penalised.
- **Symmetries.** A symmetric assembly is built from a single subunit, and only the parameters of the symmetry are refined.
- **The optimization.** Each step moves a body or a symmetry at random, and the move is kept if it lowers the χ².

The whole procedure is described by a short configuration script, so you decide which parts move, how far, and for how long.

> [!NOTE]
> Rigid-body refinement is still in active development.

## Example
Here we build the dimer of the SARS-CoV-2 main protease from a single chain, [SASDJG5_single.pdb](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/SASDJG5_single.pdb), by refining a two-fold symmetry against a SAXS measurement of the dimer, [SASDJG5.dat](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/SASDJG5.dat). Save the following script as `dimer.conf`. It declares the symmetry, keeps the two copies within 30 Å of each other, and runs 1000 steps:

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

The refined dimer is written to `output/dimer/dimer.pdb`. It fits the data better than the deposited dimer model does.

![Before refinement, the two copies of the monomer sit too far apart and fit the data poorly; the refined dimer follows the measured data across the full q range.](assets/figures/rigidbody-fit.svg#only-light){ width="100%" }
![Before refinement, the two copies of the monomer sit too far apart and fit the data poorly; the refined dimer follows the measured data across the full q range.](assets/figures/rigidbody-fit-dark.svg#only-dark){ width="100%" }

For a walkthrough in the graphical interface, see the [rigid-body tutorial](https://sastutorials.org/RigidBody/RigidBodyGUI.html). The full documentation for the `rigidbody` program, including every element of the configuration scripts, is available [here](rigidbody.md).
