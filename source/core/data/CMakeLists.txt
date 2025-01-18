target_sources(ausaxs_core PRIVATE 
	"Body.cpp"
	"Molecule.cpp"

	"symmetry/BodySymmetryFacade.cpp"
	"symmetry/SymmetryManagerMT.cpp"
	"symmetry/PredefinedSymmetries.cpp"

	"state/BoundSignaller.cpp"
	"state/StateManager.cpp"
)
