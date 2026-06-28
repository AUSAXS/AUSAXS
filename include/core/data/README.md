This folder contains the data structures representing the physical system being modelled — atoms, bodies, and the molecules built from them.

Contents
- `Molecule.h` — top-level container holding one or more bodies plus their hydration shell. The primary object passed around the calculation pipeline.
- `Body.h` — the smallest collection of atoms AUSAXS works with. A body is simply one convenient "part" of a molecule (for instance a single domain); a molecule may consist of several. There is no physical criterion for the split — it is purely an organizational abstraction.
- `atoms/` — the atom types: bare `Atom`, form-factor-aware `AtomFF`, and `Water`; plus `AtomMetaData`, optional per-atom metadata (e.g. C-alpha backbone classification, occupancy) attached to a `Body` when enabled.
- `state/` — the state-manager and signaller machinery that tracks which bodies have changed, allowing histograms to be recomputed incrementally.
- `symmetry/` — symmetry descriptions (cyclic, point, predefined) and the facades that apply them to bodies and molecules.
