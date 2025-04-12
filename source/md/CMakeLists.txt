add_library(ausaxs_md OBJECT 
	"programs/gmx.cpp"
	"programs/mdrun.cpp"
	"programs/saxsmdrun.cpp"

	"programs/options/forcefields/IForcefield.cpp"
	"programs/options/water_models/IWaterModel.cpp"
	
	"simulate/FrameAnalysis.cpp"
	"simulate/SimulateBuffer.cpp"
	"simulate/SimulateMolecule.cpp"
	"simulate/SimulateSAXS.cpp"
	"simulate/TimeAnalysis.cpp"
	
	"utility/Protein.cpp"
	"utility/files/GROFile.cpp"
	"utility/files/ITPFile.cpp"
	"utility/files/MDPCreator.cpp"
	"utility/files/NDXFile.cpp"
	"utility/files/TOPFile.cpp"
	"utility/files/TopologyFile.cpp"
)

target_include_directories(ausaxs_md PUBLIC "${CMAKE_SOURCE_DIR}/include/md")
target_link_libraries(ausaxs_md ausaxs_core)
set_target_properties(ausaxs_md PROPERTIES POSITION_INDEPENDENT_CODE 1)
