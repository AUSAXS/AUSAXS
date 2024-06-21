add_library(ausaxs_em OBJECT 
	"Box.cpp"
	"Image.cpp"
	"ImageStack.cpp"
	"ObjectBounds2D.cpp"
	"ObjectBounds3D.cpp"

	"detail/ImageStackBase.cpp"
	"detail/EMFit.cpp"
	"detail/header/DummyHeader.cpp"
	"detail/header/HeaderFactory.cpp"
	"detail/header/MapHeader.cpp"
	"detail/header/MRCHeader.cpp"
	"detail/header/RECHeader.cpp"
	"detail/header/data/HeaderData.cpp"
	"detail/header/data/MRCData.cpp"
	"detail/header/data/RECData.cpp"
	
	"manager/ProteinManager.cpp"
	"manager/ProteinManagerFactory.cpp"
	"manager/SimpleProteinManager.cpp"
	"manager/SmartProteinManager.cpp"
	
	"plots/PlotImage.cpp"
)

target_include_directories(ausaxs_em PUBLIC "${CMAKE_SOURCE_DIR}/include/em")
target_link_libraries(ausaxs_em ausaxs_core)
set_target_properties(ausaxs_em PROPERTIES POSITION_INDEPENDENT_CODE 1)
