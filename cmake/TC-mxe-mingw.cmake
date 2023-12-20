set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_C_COMPILER   /home/krelle/projects/mxe/usr/bin/x86_64-w64-mingw32.static-gcc)
set(CMAKE_CXX_COMPILER /home/krelle/projects/mxe/usr/bin/x86_64-w64-mingw32.static-g++)

# where is the target environment located
set(CMAKE_FIND_ROOT_PATH  /home/krelle/projects/mxe/usr/x86_64-w64-mingw32/x86_64-w64-mingw32.static)

# adjust the default behavior of the FIND_XXX() commands:
# search programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)

# search headers and libraries in the target environment
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
