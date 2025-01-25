target_sources(ausaxs_core PRIVATE 
	"ExplicitHydration.cpp"
	"Hydration.cpp"

	"culling/BodyCounterCulling.cpp"
	"culling/CounterCulling.cpp"
	"culling/CullingFactory.cpp"
	"culling/CullingStrategy.cpp"
	"culling/NoCulling.cpp"
	"culling/OutlierCulling.cpp"
	"culling/RandomCulling.cpp"
	
	"generation/AxesHydration.cpp"
	"generation/GridBasedHydration.cpp"
	"generation/HydrationFactory.cpp"
	"generation/JanHydration.cpp"
	"generation/NoHydration.cpp"
	"generation/PepsiHydration.cpp"
	"generation/RadialHydration.cpp"
)
