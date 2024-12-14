target_sources(ausaxs_core PRIVATE 
	"Body.cpp"
	"Molecule.cpp"
	"SymmetryManager.cpp"

	"detail/AtomCollection.cpp"

	"record/PDBAtom.cpp"
	"record/Footer.cpp"
	"record/Header.cpp"
	"record/Record.cpp"
	"record/Terminate.cpp"
	"record/PDBWater.cpp"

	"state/BoundSignaller.cpp"
	"state/StateManager.cpp"
	"state/UnboundSignaller.cpp"
)
