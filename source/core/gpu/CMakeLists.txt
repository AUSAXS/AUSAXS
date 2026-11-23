target_sources(ausaxs_core PRIVATE
	"GPUInstaller.cpp"
)
if (EMSCRIPTEN)
	target_sources(ausaxs_core PRIVATE "GPULoaderUnavailable.cpp")
else()
	target_sources(ausaxs_core PRIVATE "GPULoader.cpp")
endif()
