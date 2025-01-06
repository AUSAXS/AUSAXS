target_sources(ausaxs_core PRIVATE 
	"detail/CIFReader.cpp"
	"detail/PDBReader.cpp"
	"detail/PDBWriter.cpp"
	"detail/XYZWriter.cpp"

	"pdb/Footer.cpp"
	"pdb/Header.cpp"
	"pdb/PDBAtom.cpp"
	"pdb/PDBStructure.cpp"
	"pdb/PDBWater.cpp"
	"pdb/Record.cpp"
	"pdb/Terminate.cpp"

	"File.cpp"
	"Folder.cpp"
	"Reader.cpp"
	"Writer.cpp"
)
