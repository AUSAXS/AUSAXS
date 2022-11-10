add_executable(new_hydration new_hydration.cpp ${BACKWARD_ENABLE})
target_link_libraries(new_hydration saxs)

add_executable(hist hist.cpp ${BACKWARD_ENABLE})
target_link_libraries(hist saxs)

add_executable(main main.cpp ${BACKWARD_ENABLE})
target_link_libraries(main saxs)

add_executable(em em.cpp ${BACKWARD_ENABLE})
target_link_libraries(em saxs)

add_executable(optimize_radius optimize_radius.cpp ${BACKWARD_ENABLE})
target_link_libraries(optimize_radius saxs)

add_executable(intensity_fitter intensity_fitter.cpp ${BACKWARD_ENABLE})
target_link_libraries(intensity_fitter saxs)

add_executable(rigidbody rigidbody.cpp ${BACKWARD_ENABLE})
target_link_libraries(rigidbody saxs)

add_executable(stuff stuff.cpp ${BACKWARD_ENABLE})
target_link_libraries(stuff saxs)

add_executable(rotate_map rotate_map.cpp ${BACKWARD_ENABLE})
target_link_libraries(rotate_map saxs)

add_executable(consistency consistency.cpp ${BACKWARD_ENABLE})
target_link_libraries(consistency saxs)

add_executable(fit_consistency fit_consistency.cpp ${BACKWARD_ENABLE})
target_link_libraries(fit_consistency saxs)

add_executable(fit_consistency2 fit_consistency2.cpp ${BACKWARD_ENABLE})
target_link_libraries(fit_consistency2 saxs)

add_executable(rebin rebin.cpp ${BACKWARD_ENABLE})
target_link_libraries(rebin saxs)

add_executable(unit_cell unit_cell.cpp ${BACKWARD_ENABLE})
target_link_libraries(unit_cell saxs)

add_executable(em_fitter em_fitter.cpp ${BACKWARD_ENABLE})
target_link_libraries(em_fitter saxs)

add_executable(em_pdb_fitter em_pdb_fitter.cpp ${BACKWARD_ENABLE})
target_link_libraries(em_pdb_fitter saxs)

add_executable(plot_data plot_data.cpp ${BACKWARD_ENABLE})
target_link_libraries(plot_data saxs)

add_executable(inspect_map inspect_map.cpp ${BACKWARD_ENABLE})
target_link_libraries(inspect_map saxs)

add_executable(plot plot.cpp ${BACKWARD_ENABLE})
target_link_libraries(plot saxs)

if (DEBUG)
	add_backward(new_hydration)
	add_backward(hist)
	add_backward(main)
	add_backward(em)
	add_backward(optimize_radius)
	add_backward(intensity_fitter)
	add_backward(rigidbody)
	add_backward(stuff)
	add_backward(rotate_map)
	add_backward(consistency)
	add_backward(fit_consistency)
	add_backward(fit_consistency2)
	add_backward(rebin)
	add_backward(unit_cell)
	add_backward(em_fitter)
	add_backward(em_pdb_fitter)
	add_backward(plot_data)
	add_backward(inspect_map)
	add_backward(plot)
endif()
