This folder contains the voxel grid used to represent the volume occupied by a molecule.

The grid discretizes space into cells, providing a fast spatial lookup for hydration-shell placement and for estimating the excluded volume displaced by the solute.

Contents
- `Grid.h` — the main grid class: voxel storage, atom insertion, and volume queries.
- `expansion/` — strategies for expanding atoms into the grid, e.g. spherical expansion of an atom into the surrounding voxels.
- `exv/` — excluded-volume strategies that derive the displaced solvent volume from the occupied grid cells.
- `detail/` — internal helpers supporting the above.
