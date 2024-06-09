file(GLOB_RECURSE INCLUDE "${CMAKE_SOURCE_DIR}/include/md/*.h")
file(GLOB_RECURSE SOURCE "${CMAKE_SOURCE_DIR}/source/md/*.cpp")
add_library(ausaxs_gmx ${INCLUDE} ${SOURCE})
target_include_directories(gmx PUBLIC "${CMAKE_SOURCE_DIR}/include")
