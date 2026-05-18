This folder contains the rigidbody optimizer, which refines a structure by moving its bodies as rigid units.

A molecule is first divided into bodies (`BodySplitter.h`); each optimization step then selects a body, applies a transformation, regenerates the hydration shell, and accepts or rejects the move based on the fit to the SAXS data.

Contents
- `constraints/` — constraints linking bodies together (distance, overlap, attractors, repellers) and the manager enforcing them.
- `selection/` — strategies for choosing which body or constraint to perturb next.
- `parameters/` — generation of the transformation parameters for each step, including decay of the step size over time.
- `transform/` — the transformations themselves and how they propagate through linked bodies.
- `controller/` — the optimization loop and acceptance criterion (e.g. Metropolis).
- `sequencer/` — a high-level scripting interface for assembling a custom optimization run from reusable elements.
