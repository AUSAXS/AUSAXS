---
title: Compiling AUSAXS from source - AUSAXS
description: Compile the AUSAXS small-angle scattering calculator from source with CMake on Linux, macOS or Windows.
---

# Compiling from source

The [Python package](installation.md) is the recommended way to install AUSAXS, but compiling it yourself is also fairly straightforward. It only requires the `curl` library and its development headers. Since the codebase is written in modern `C++20`, a recent compiler supporting this is also required. The compilation has been verified to work with GCC v11+, Clang v15+, and MSVC 2022+.

## Linux
1. Make sure you have the prerequisites installed  
`sudo apt-get install cmake make g++ git libcurl4-openssl-dev`

2. Clone the repository  
`git clone https://github.com/AUSAXS/AUSAXS.git`

3. Configure the build from the root directory  
`cmake -B build -S . -DCMAKE_BUILD_TYPE=Release`

4. Compile the executable  
`cmake --build build --target ausaxs -jX`  
where X is the number of worker threads. The executable is placed in `build/bin/`. 

5. (optional) To convert the output `.plot` files to images, you need a few Python packages.  
`python -m pip install matplotlib numpy scipy`  

## macOS
1. Make sure both CMake and a recent compiler, e.g. Xcode or the Xcode Command Line Tools, are available on your system  
`brew install cmake`

2. Follow steps 2 through 5 from the Linux guide. 

## Windows
1. Make sure CURL is available on your system, e.g. through vcpkg  
`vcpkg install curl --triplet x64-windows-static`

2. Download or clone the repository  
`git clone https://github.com/AUSAXS/AUSAXS.git`

3. Configure the build with vcpkg's toolchain file, so CMake can find CURL  
`cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=<vcpkg root>/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows-static`

4. Compile the executable  
`cmake --build build --target ausaxs --config Release -jX`  
The executable is placed in `build/bin/Release/`. Note that this can be very memory-intensive with the MSVC compiler. 

5. (optional) To convert the output `.plot` files to images, you will also need `Python` and a few additional packages.  
`python -m pip install matplotlib numpy scipy`

## Dependencies
The following are automatically fetched by CMake:

* [Generalized Constant Expression Math](https://github.com/kthohr/gcem), for constexpr math functions on compilers other than GCC. `-DUSE_SYSTEM_GCEM=ON`
* [BSThreadPool](https://github.com/bshoshany/thread-pool), for multithreading. `-DUSE_SYSTEM_THREADPOOL=ON`
* [CLI11](https://github.com/CLIUtils/CLI11), for interpreting command line arguments. `-DUSE_SYSTEM_CLI11=ON`
* [pocketfft](https://github.com/mreineck/pocketfft), for the Fourier transform of the grid-based excluded volume. `-DUSE_SYSTEM_POCKETFFT=ON`
* [miniz](https://github.com/richgel999/miniz), for unpacking the GPU driver in `ausaxs gpu install`. `-DUSE_SYSTEM_MINIZ=ON`
* [Catch2](https://github.com/catchorg/Catch2), for running tests. `-DUSE_SYSTEM_CATCH=ON`, or `-DBUILD_TESTS=OFF` to skip the tests entirely.

Each of these can be replaced by a system-provided version with the CMake option listed next to it.

## Using your build from Python
The Python package bundles its own copy of the AUSAXS library. To make it use one you compiled yourself instead, e.g. to try your changes to the library in the [graphical interface](../documentation/cli/fit.md#graphical-user-interface), build the shared library:
```bash
cmake --build build --target libausaxs -jX
```
and point the package at it with `ausaxs setup`:
```bash
ausaxs setup --relink build/lib/libausaxs.so
```
The library is named `libausaxs.dylib` on macOS, and on Windows it is `libausaxs.dll` in `build/bin/Release/`. The path is checked by loading the library, and remembered across sessions, so rebuilding the library is enough to pick up later changes. The other options are:

| Command | Effect |
|---|---|
| `ausaxs setup` | Show which library is in use, and where it comes from. |
| `ausaxs setup --reset` | Go back to the bundled library. |

For a one-off run, the `AUSAXS_LIB` environment variable can be set to the library path instead. A relinked path takes precedence over it.

## Next steps
Check that the build works by running `build/bin/ausaxs --version`, then continue with the [next steps](installation.md#next-steps).
