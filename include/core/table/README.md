This folder contains lookup tables that cache expensive repeated computations.

The principal one is the Debye table, which stores the `sin(qd)/(qd)` factor for every combination of scattering vector `q` and distance bin `d`. Pre-computing this table turns the histogram-to-intensity conversion into a simple weighted sum. Both array-backed and vector-backed implementations are provided, selected via `DebyeTableManager.h`.
