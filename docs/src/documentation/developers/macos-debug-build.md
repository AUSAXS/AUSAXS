---
title: macOS debug build - AUSAXS
description: How to compile AUSAXS in debug mode on macOS.
---

# macOS debug build

Ensure `libiberty` is included in the binutils brew by adding `--enable-install-libiberty` to its configuration:  
`brew edit binutils`

Then rebuild it from scratch:  
`brew reinstall -s binutils`

The Xcode distribution of binutils will take precedence in the CMake configuration unless we specifically tell it otherwise.  
`export LDFLAGS="-L/usr/local/opt/binutils/lib"`  
`export CPPFLAGS="-I/usr/local/opt/binutils/include"`

Finally, CMake can be invoked in debug mode as:  
`cmake -B build -S . -DCMAKE_BUILD_TYPE=debug`

This should allow all targets to be built in debug mode. 
