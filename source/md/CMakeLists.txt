add_library(ausaxs_md OBJECT 
	"programs/gmx.cpp"
	"programs/mdrun.cpp"
	"programs/saxsmdrun.cpp"
	
	"simulate/buffer.cpp"
	"simulate/frameanalysis.cpp"
	"simulate/molecule.cpp"
	"simulate/saxs.cpp"
	"simulate/timeanalysis.cpp"
	
	"utility/Protein.cpp"
	"utility/files/GROFile.cpp"
	"utility/files/ITPFile.cpp"
	"utility/files/MDPCreator.cpp"
	"utility/files/TOPFile.cpp"
	"utility/files/TopologyFile.cpp"
)

target_include_directories(ausaxs_md PUBLIC "${CMAKE_SOURCE_DIR}/include/md")
target_link_libraries(ausaxs_md ausaxs_core)
set_target_properties(ausaxs_md PROPERTIES POSITION_INDEPENDENT_CODE 1)
