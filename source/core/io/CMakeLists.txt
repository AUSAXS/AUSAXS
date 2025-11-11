target_sources(ausaxs_core PRIVATE 
	"detail/structure/CIFReader.cpp"
	"detail/structure/PDBReader.cpp"
	"detail/structure/XYZReader.cpp"
	"detail/structure/PDBWriter.cpp"
	"detail/trajectory/XYZWriter.cpp"

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
