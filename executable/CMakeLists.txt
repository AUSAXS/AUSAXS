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

add_executable(rigidbody rigidbody.cpp)
target_link_libraries(rigidbody saxs)

add_executable(stuff stuff.cpp)
target_link_libraries(stuff saxs)

add_executable(rotate_map rotate_map.cpp)
target_link_libraries(rotate_map saxs)

add_executable(consistency consistency.cpp)
target_link_libraries(consistency saxs)

add_executable(fit_consistency fit_consistency.cpp)
target_link_libraries(fit_consistency saxs)

add_executable(rebin rebin.cpp)
target_link_libraries(rebin saxs)

add_executable(unit_cell unit_cell.cpp)
target_link_libraries(unit_cell saxs)
