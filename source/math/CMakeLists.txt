add_library(ausaxs_math OBJECT 
	"CubicSpline.cpp"
	"detail/Diagnostics.cpp"
	"detail/Format.cpp"
	"detail/Streams.cpp"
	"LUPDecomposition.cpp"
	"MatrixUtils.cpp"
	"MovingAverager.cpp"
	"PeakFinder.cpp"
	"Superposition.cpp"
	"SymmetricEigen.cpp"
)

target_include_directories(ausaxs_math PUBLIC "${CMAKE_SOURCE_DIR}/include/math")
set_target_properties(ausaxs_math PROPERTIES POSITION_INDEPENDENT_CODE 1)
