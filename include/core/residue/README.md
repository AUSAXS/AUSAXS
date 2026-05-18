This folder contains the residue storage used to assign atomic properties.

`ResidueStorage.h` holds per-residue information — such as the number of bound hydrogens for each atom — which is needed to determine effective form factors. Definitions are looked up by residue name and cached, with missing residues fetched on demand.
