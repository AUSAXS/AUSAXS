This folder contains the standalone executables built for AUSAXS releases.

Binaries
- `saxs_fitter` — SAXS fitting command-line tool (see `executable/saxs_fitter.cpp`).
- `em_fitter` — electron-microscopy fitting utility (see `executable/em_fitter.cpp`).
- `rigidbody_optimizer` — rigid-body configuration optimizer (see `executable/rigidbody_optimizer.cpp`).
- `saxs_fitter_gui` — GUI version of the SAXS fitting tool (see `executable/gui/saxs_fitter_gui.cpp`)
- `em_fitter_gui` — GUI version of the electron-microscopy fitting tool (see `executable/gui/em_fitter_gui.cpp`)

Build
- All executables can be built by filename. 
  1. `cmake -B build -S .`
  2. `cmake --build build --target <filename> -jX`
- This will build the binaries and place them in the build output directory (`build/bin`).

Usage
- Examples and CLI option documentation are maintained on the project wiki: https://github.com/AUSAXS/AUSAXS/wiki

Support
- Report issues: https://github.com/AUSAXS/AUSAXS/issues
- Contributing: see [`CONTRIBUTING.md`](/CONTRIBUTING.md)
- License: see [`LICENSE`](/LICENSE)