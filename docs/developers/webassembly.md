# Pyodide backend

The manually triggered **Build Pyodide backend** workflow builds `libausaxs.so`, an Emscripten side module intended to be bundled with the pyAUSAXS Python wheel. It targets the [pyemscripten_2026_0 ABI](https://pyodide.org/en/stable/development/abi/314.html), used by Pyodide with Python 3.14: Emscripten 5.0.3 and WebAssembly exception handling. It is not a native Linux library or a standalone JavaScript module.

To build locally, activate Emscripten 5.0.3 and run:

```sh
emcmake cmake -S . -B build-wasm -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTS=OFF -DBUILD_EXECUTABLES=OFF -DGPU_DOWNLOADER=OFF
cmake --build build-wasm --target libausaxs -j2
```

The resulting `build-wasm/lib/libausaxs.so` is single-threaded, without cURL, a GPU backend, or shell commands. PocketFFT is enabled by default and can be disabled with `-DPOCKETFFT=OFF`. Built-in residue definitions are available, but unrecognized residues cannot be downloaded; provide any additional required residue files in the virtual cache. Calls execute synchronously, so a Pyodide Web Worker is appropriate for long calculations.

## Integrating with pyAUSAXS

The workflow artifact is the backend only; packaging it into a wheel and loading it in the wiki requires corresponding changes in pyAUSAXS:

1. Put `libausaxs.so` in `pyausaxs/resources/` when building a Pyodide-targeted wheel. Tag the wheel for the matching Pyodide platform instead of publishing it as a pure Python wheel. The wiki must load a Pyodide release using the same ABI and install that wheel and its dependencies.
2. Teach the Python loader that `sys.platform == "emscripten"` uses the bundled `.so` (the current OS detection treats it as unknown). Load it with `ctypes.CDLL` from the Pyodide filesystem.
3. Skip the native CPU feature check and replace the multiprocessing-based integration check with an in-process `test_integration` call on Pyodide, where subprocesses are unavailable.

The pyAUSAXS wheel build, loader, and dependency compatibility are maintained separately; this repository's workflow does not modify or publish the wheel. Other Pyodide platforms require their own matching Emscripten version and exception ABI.
