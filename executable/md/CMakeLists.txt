add_executable(simulate simulate.cpp)
target_link_libraries(simulate ausaxs_core ausaxs_math ausaxs_md CLI11::CLI11)

add_executable(simulate_buffer simulate_buffer.cpp)
target_link_libraries(simulate_buffer ausaxs_core ausaxs_math ausaxs_md CLI11::CLI11)

add_executable(time_analysis time_analysis.cpp)
target_link_libraries(time_analysis ausaxs_core ausaxs_math ausaxs_md CLI11::CLI11)

add_executable(frame_analysis frame_analysis.cpp)
target_link_libraries(frame_analysis ausaxs_core ausaxs_math ausaxs_md CLI11::CLI11)

add_executable(extract_profile_water extract_profile_water.cpp)
target_link_libraries(extract_profile_water ausaxs_core ausaxs_math ausaxs_md CLI11::CLI11)

add_executable(background_schemes background_schemes.cpp)
target_link_libraries(background_schemes ausaxs_core ausaxs_math ausaxs_md CLI11::CLI11)