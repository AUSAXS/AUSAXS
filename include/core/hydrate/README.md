This folder contains the hydration-shell models — the layer of ordered water surrounding the solute that contributes significantly to the SAXS signal.

Contents
- `Hydration.h` and the `ExplicitHydration`, `ImplicitHydration`, and `EmptyHydration` variants — the different ways a hydration shell can be represented, from individually placed water molecules to no shell at all.
- `generation/` — strategies that place water molecules around the molecule, typically using the voxel grid.
- `culling/` — strategies that thin out an over-dense generated shell down to a physically reasonable coverage.
