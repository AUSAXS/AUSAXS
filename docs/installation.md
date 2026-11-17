---
title: Installing AUSAXS
description: Install the AUSAXS small-angle scattering calculator with pip install pyausaxs, download a precompiled binary, or compile it from source on Linux, macOS or Windows.
---

# Installing AUSAXS

`AUSAXS` can be installed in three ways:  
1. [Installing the Python wrapper](#installation-via-the-python-wrapper),  
2. [Downloading the latest release](#downloading-the-latest-github-release), and  
3. [Compiling it yourself](#compiling-the-library-yourself). 

## Installation via the Python wrapper

The latest version can be easily installed by running the following: 
```python
pip install pyausaxs[all]
```

See the [pyausaxs repository](https://github.com/AUSAXS/pyAUSAXS) for more information. 

## Downloading the latest GitHub release

Recent precompiled executables for all platforms are available in the [releases](https://github.com/AUSAXS/AUSAXS/releases) section.

## Compiling the library yourself

If you prefer compiling the library yourself, this is also fairly straightforward, requiring only that `curl` is available on the system. Since the codebase is written in modern `C++20`, a recent compiler supporting this is also required. The compilation has been verified to work with GCC v11+, Clang v15+, and MSVC 2022+.

### Linux
1. Make sure you have the prerequisites installed  
`apt-get install cmake make g++ libcurl4-openssl-dev`

2. Clone this repository  
`git clone https://github.com/AUSAXS/AUSAXS.git`.

3. Run the build command from the root directory  
`cmake -B build -S .`

4. Compile the executable  
`cmake --build build --target ausaxs -jX`  
where X is the number of worker threads. The executable is placed in `build/bin/`. 

5. (optional) To convert the output `.plot` files to images, you need a few Python packages.  
`pip3 install matplotlib numpy scipy`  

### MacOS
1. Make sure both CMake and a recent compiler, e.g. Xcode, is available on your system  
`brew install cmake`

2. Follow step 2 through 5 from the Linux guide. 

### Windows
1. Make sure CURL is available on your system, e.g. through vcpkg  
`vcpkg install curl`

2. Download or clone this repository  
`git clone https://github.com/AUSAXS/AUSAXS.git`.

3. Compile the `ausaxs` target as in step 3 and 4 of the Linux guide. Note that this can be very memory-intensive with the MSVC compiler. 

4. (optional) To convert the output `.plot` files to images, you will also need `Python` and a few additional packages.  
`pip3 install matplotlib numpy scipy`

### Dependencies
The following are automatically fetched by CMake:
*	[Generalized Constant Expression Math](https://github.com/kthohr/gcem), for constexpr math support beyond GCC. 
*	[BSThreadPool](https://github.com/bshoshany/thread-pool), for multithreading. 
*	[CLI11](https://github.com/CLIUtils/CLI11), for interpreting command line arguments. 
*	[pocketfft](https://github.com/mreineck/pocketfft), for the Fourier transform of the grid-based excluded volume. 
*	[Catch2](https://github.com/catchorg/Catch2), for running tests.

Each of these can be replaced by a system-provided version with the corresponding `-DUSE_SYSTEM_<NAME>=ON` CMake option. 
