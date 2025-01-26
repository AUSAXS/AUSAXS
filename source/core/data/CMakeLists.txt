target_sources(ausaxs_core PRIVATE 
	"Body.cpp"
	"Molecule.cpp"

	"symmetry/BodySymmetryFacade.cpp"
	"symmetry/MoleculeSymmetryFacade.cpp"
	"symmetry/SymmetryManagerMT.cpp"
	"symmetry/PredefinedSymmetries.cpp"
	"symmetry/detail/SymmetryHelpers.cpp"

	"state/BoundSignaller.cpp"
	"state/StateManager.cpp"
)
