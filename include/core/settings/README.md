This folder contains the global settings that configure an AUSAXS run.

Settings are grouped by area (`GridSettings`, `FitSettings`, `EMSettings`, `HistogramSettings`, and so on), each exposing a set of tunable values with sensible defaults. The remaining headers handle reading settings from a configuration file, registering them, and validating them. `All.h` is a convenience header pulling in every settings group.
