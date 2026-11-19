---
title: Compiling AUSAXS from source - AUSAXS
description: Compile the AUSAXS small-angle scattering calculator from source with CMake on Linux, macOS or Windows.
---

# Compiling from source

The [Python package](installation.md) is the recommended way to install AUSAXS, but compiling it yourself is also fairly straightforward. It only requires that `curl` is available on the system. Since the codebase is written in modern `C++20`, a recent compiler supporting this is also required. The compilation has been verified to work with GCC v11+, Clang v15+, and MSVC 2022+.

## Linux
1. Make sure you have the prerequisites installed  
`apt-get install cmake make g++ libcurl4-openssl-dev`

2. Clone the repository  
`git clone https://github.com/AUSAXS/AUSAXS.git`.

3. Run the build command from the root directory  
`cmake -B build -S .`

4. Compile the executable  
`cmake --build build --target ausaxs -jX`  
where X is the number of worker threads. The executable is placed in `build/bin/`. 

5. (optional) To convert the output `.plot` files to images, you need a few Python packages.  
`pip3 install matplotlib numpy scipy`  

## macOS
1. Make sure both CMake and a recent compiler, e.g. Xcode, is available on your system  
`brew install cmake`

2. Follow step 2 through 5 from the Linux guide. 

## Windows
1. Make sure CURL is available on your system, e.g. through vcpkg  
`vcpkg install curl`

2. Download or clone the repository  
`git clone https://github.com/AUSAXS/AUSAXS.git`.

3. Configure the build with vcpkg's toolchain file, so CMake can find CURL  
`cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=<vcpkg root>/scripts/buildsystems/vcpkg.cmake`

4. Compile the `ausaxs` target as in step 4 of the Linux guide. Note that this can be very memory-intensive with the MSVC compiler. 

5. (optional) To convert the output `.plot` files to images, you will also need `Python` and a few additional packages.  
`pip3 install matplotlib numpy scipy`

## Dependencies
The following are automatically fetched by CMake:

* [Generalized Constant Expression Math](https://github.com/kthohr/gcem), for constexpr math support beyond GCC. 
* [BSThreadPool](https://github.com/bshoshany/thread-pool), for multithreading. 
* [CLI11](https://github.com/CLIUtils/CLI11), for interpreting command line arguments. 
* [pocketfft](https://github.com/mreineck/pocketfft), for the Fourier transform of the grid-based excluded volume. 
* [Catch2](https://github.com/catchorg/Catch2), for running tests.

Each of these can be replaced by a system-provided version with the corresponding `-DUSE_SYSTEM_<NAME>=ON` CMake option.

## Using your build from Python
The Python package bundles its own copy of the AUSAXS library. To make it use one you compiled yourself instead, e.g. to try out changes in the [graphical interface](../documentation/cli/fit.md#graphical-user-interface), build the shared library:
```bash
cmake --build build --target libausaxs -jX
```
and point the package at it with `ausaxs setup`:
```bash
ausaxs setup --relink build/lib/libausaxs.so
```
The library is named `libausaxs.dylib` on macOS and `libausaxs.dll` on Windows. The path is checked by loading the library, and remembered across sessions, so rebuilding the library is enough to pick up later changes. The other options are:

| Command | Effect |
|---|---|
| `ausaxs setup` | Show which library is in use, and where it comes from. |
| `ausaxs setup --reset` | Go back to the bundled library. |

For a one-off run, the `AUSAXS_LIB` environment variable can be set to the library path instead. A relinked path takes precedence over it.

## Next steps
Check that the build works by running `build/bin/ausaxs --version`, then continue with the [next steps](installation.md#next-steps).
