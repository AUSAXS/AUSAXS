target_sources(ausaxs_core PRIVATE 
	"Body.cpp"
	"Molecule.cpp"

	"detail/AtomCollection.cpp"

	"record/Atom.cpp"
	"record/Footer.cpp"
	"record/Header.cpp"
	"record/Record.cpp"
	"record/Terminate.cpp"
	"record/Water.cpp"

	"state/BoundSignaller.cpp"
	"state/StateManager.cpp"
	"state/UnboundSignaller.cpp"
)
