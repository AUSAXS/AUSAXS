add_library(ausaxs_crystal OBJECT 
	"CrystalScattering.cpp"
	"Fval.cpp"

	"io/CrystalReaderFactory.cpp"
	"io/GridReader.cpp"
	"io/PDBReader.cpp"
	"io/UnitCellReader.cpp"

	"miller/AllMillers.cpp"
	"miller/FibonacciMillers.cpp"
	"miller/Miller.cpp"
	"miller/MillerGenerationFactory.cpp"
	"miller/ReducedMillers.cpp"
)

target_include_directories(ausaxs_crystal PUBLIC "${CMAKE_SOURCE_DIR}/include/crystal")
target_link_libraries(ausaxs_crystal PRIVATE ausaxs_core ausaxs_math)
set_target_properties(ausaxs_crystal PROPERTIES
	PROPERTIES_OUTPUT_NAME ausaxs_crystal
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)
