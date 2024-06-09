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

set_target_properties(ausaxs_crystal PROPERTIES
	PROPERTIES_OUTPUT_NAME ausaxs_crystal
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)
