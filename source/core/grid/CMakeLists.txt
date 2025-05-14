target_sources(ausaxs_core PRIVATE 
	"Grid.cpp"
	
	"detail/GridMember.cpp"
	"detail/GridObj.cpp"
	"detail/GridSurfaceDetection.cpp"
	"detail/RadialLineGenerator.cpp"

	"exv/GridExvStrategy.cpp"
	"exv/RawGridExv.cpp"
	"exv/RawGridWithSurfaceExv.cpp"
)