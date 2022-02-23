add_executable(new_hydration new_hydration.cpp)
target_link_libraries(new_hydration saxs)

add_executable(hist hist.cpp)
target_link_libraries(hist saxs)

add_executable(main main.cpp)
target_link_libraries(main saxs)

add_executable(em em.cpp)
target_link_libraries(em saxs)

add_executable(optimize_radius optimize_radius.cpp)
target_link_libraries(optimize_radius saxs)

add_executable(intensity_fitter intensity_fitter.cpp)
target_link_libraries(intensity_fitter saxs)

add_executable(rigidbody_opt rigidbody.cpp)
target_link_libraries(rigidbody_opt saxs)
