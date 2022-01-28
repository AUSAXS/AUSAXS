#include "io/PDBWriter.h"
#include "io/Writer.h"
#include "io/File.h"
#include "data/Terminate.h"
#include "data/Atom.h"
#include "data/Hetatom.h"

#include <fstream>
#include <algorithm>

void PDBWriter::write(const string& output_path) {
    file->refresh();
    std::ofstream output(output_path);
    if (!output.is_open()) {throw std::ios_base::failure("Error in PDB_file::write: Could not open file \"" + output_path + "\"");}
    output << as_pdb() << std::flush;
    output.close();
    cout << "Output written to file " + output_path + "." << endl;
}

string PDBWriter::as_pdb() const {
    File& f = *file;
    string s;
    s += f.header.get();

    unsigned int i_ter = f.terminate.serial;
    bool printed_ter = false;
    for (unsigned int i = 0; i < f.protein_atoms.size(); i++) {
        if (i == i_ter) { // check if this is where the terminate is supposed to go
            s += f.terminate.as_pdb(); // write it if so
            printed_ter = true;
        }
        s += f.protein_atoms[i].as_pdb();
    }

    // print terminate if missing
    if (!printed_ter) {s += f.terminate.as_pdb();}

    // print hetatoms
    std::for_each(f.hydration_atoms.begin(), f.hydration_atoms.end(), [&s] (const Hetatom& atom) {s += atom.as_pdb();});

    s += f.footer.get();
    return s;
}