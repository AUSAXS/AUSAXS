# Building for WebAssembly

The Emscripten build produces a self-contained JavaScript loader and WebAssembly module exposing the AUSAXS C API. It runs synchronously on one thread, without cURL, a GPU backend, or command-line tools. It does not require browser cross-origin isolation.

## Build and test

Install and activate the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html). This build has been tested with Emscripten 7.0.1 and Node.js 24.

```sh
emcmake cmake -S . -B build-wasm
cmake --build build-wasm -j
ctest --test-dir build-wasm --output-on-failure
```

The output is `build-wasm/bin/ausaxs.js` and `build-wasm/bin/ausaxs.wasm`. The tests use Node, which `emcmake` normally sets as the CMake cross-compiling emulator. To omit tests, configure with `-DBUILD_TESTS=OFF`. PocketFFT is enabled by default and can be disabled with `-DPOCKETFFT=OFF`.

## Use in a browser

Serve both files from the same directory over HTTP(S), with the `.wasm` file served as `application/wasm`. Load `ausaxs.js` as a classic script and await its factory:

```html
<script src="ausaxs.js"></script>
<script>
async function example() {
    const ausaxs = await createAUSAXS();
    const value = ausaxs._malloc(4);
    try {
        ausaxs.HEAP32[value / 4] = 41;
        ausaxs._test_integration(value);
        console.log(ausaxs.HEAP32[value / 4]); // 42
    } finally {
        ausaxs._free(value);
    }
}
example().catch(console.error);
</script>
```

The [C API headers](../../include/api/api/) describe the exported functions and their pointer, status, and object-lifetime conventions; CLI entry points are not included. JavaScript can call `_function_name` directly, or use `ccall` / `cwrap`. `_malloc`, `_free`, `HEAP32`, `HEAPU32`, `HEAPF64`, `UTF8ToString`, `addFunction`, `removeFunction`, and `FS` are available on the module. WebAssembly memory can grow, so obtain a fresh `HEAP*` view after a call that may allocate. The [smoke test](../../tests/wasm/smoke.js) demonstrates numerical calls, callbacks, status checking, object deallocation, error messages, and loading a PDB into `FS`.

Files must be supplied through Emscripten's virtual filesystem, for example `ausaxs.FS.writeFile("/input.pdb", bytes)`. Its default in-memory filesystem is not persistent across page reloads. Downloads are disabled: built-in residue definitions are available, but unrecognized residues cannot be downloaded and can affect results. Set up any additional required residue files in the virtual cache beforehand.

The `threads` setting is fixed at one and GPU use is disabled; attempts to enable them are rejected with a warning. Shell commands are unavailable. Calculations execute on the calling JavaScript thread, so use small examples or run the module in a Web Worker to keep the page responsive.

## Python

This output is an Emscripten JavaScript/WebAssembly module, **not** a Python runtime or a Python wheel. Using the separate pyAUSAXS package inside a browser Python runtime requires additional integration: its desktop `ctypes` loader and subprocess-based integration check cannot simply load this module as a native shared library. Choose a browser Python runtime and adapt its bridge and packaging separately.
