add_executable(new_hydration new_hydration.cpp ${sources})
add_executable(hist hist.cpp ${sources})
add_executable(main main.cpp ${sources})
add_executable(optimize_radius optimize_radius.cpp ${sources})
add_executable(intensity_fitter intensity_fitter.cpp ${sources})

target_link_libraries(new_hydration ${Boost_LIBRARIES})
target_link_libraries(hist ${Boost_LIBRARIES})
