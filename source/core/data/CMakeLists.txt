target_sources(ausaxs_core PRIVATE 
	"Body.cpp"
	"Molecule.cpp"

	"symmetry/BodySymmetryFacade.cpp"
	"symmetry/MoleculeSymmetryFacade.cpp"
	"symmetry/PointSymmetry.cpp"
	"symmetry/PredefinedSymmetries.cpp"
	"symmetry/CyclicSymmetry.cpp"

	"state/BoundSignaller.cpp"
	"state/StateManager.cpp"
)
