target_sources(ausaxs_core PRIVATE 
	"Axis.cpp"
	"Axis3D.cpp"
	"Console.cpp"
	"Exceptions.cpp"
	"Limit.cpp"
	"Limit2D.cpp"
	"Limit3D.cpp"
	"Logging.cpp"
	"StringUtils.cpp"
	"Type.cpp"
	"Utility.cpp"
)

if (EMSCRIPTEN)
	target_sources(ausaxs_core PRIVATE "CurlUnavailable.cpp" "MultiThreadingUnavailable.cpp")
else()
	target_sources(ausaxs_core PRIVATE "Curl.cpp" "MultiThreading.cpp")
endif()
