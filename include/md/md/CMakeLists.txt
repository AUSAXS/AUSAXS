cmake_minimum_required(VERSION 3.15)

file(GLOB_RECURSE INCLUDE "${CMAKE_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE SOURCE "${CMAKE_SOURCE_DIR}/source/*.cpp")
add_library(gmx ${INCLUDE} ${SOURCE})
target_include_directories(gmx PUBLIC "${CMAKE_SOURCE_DIR}/include")
#target_link_libraries(gmx CLI11::CLI11 -static-libgcc -static-libstdc++)
target_link_libraries(gmx CLI11::CLI11)


