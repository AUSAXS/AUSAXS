target_sources(ausaxs_core PRIVATE 
	"Body.cpp"
	"Molecule.cpp"

	"symmetry/BodySymmetryFacade.cpp"
	"symmetry/MoleculeSymmetryFacade.cpp"
	"symmetry/ISymmetry.cpp"
	"symmetry/PointSymmetry.cpp"
	"symmetry/PredefinedSymmetries.cpp"
	"symmetry/CyclicSymmetry.cpp"
	"symmetry/IcosahedralSymmetry.cpp"
	"symmetry/IPolyhedralSymmetry.cpp"
	"symmetry/OctahedralSymmetry.cpp"
	"symmetry/CompositeSymmetry.cpp"
	"symmetry/ReferenceSymmetry.cpp"
	"symmetry/PairScheduleUtil.cpp"
	"symmetry/TetrahedralSymmetry.cpp"
	"symmetry/DihedralSymmetry.cpp"

	"state/BoundSignaller.cpp"
	"state/StateManager.cpp"
)
