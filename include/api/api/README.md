This folder contains the external API exposing AUSAXS functionality to other languages and tools.

Contents
- `pyausaxs/` — the C interface backing the [pyausaxs](https://github.com/AUSAXS/pyAUSAXS) Python wrapper, covering molecules, fitting, form factors, I/O, and settings.
- `cli/` — entry points used by the command-line executables.
- `api_sasview.h` — integration hooks for SasView.
- `ObjectStorage.h` — keeps C++ objects alive across the language boundary, handing out opaque handles to callers.
