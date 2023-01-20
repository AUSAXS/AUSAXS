#include <io/PDBWriter.h>
#include <io/Writer.h>
#include <io/File.h>
#include <data/Terminate.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

#include <fstream>
#include <algorithm>

void PDBWriter::write(std::string output_path) {
    file->refresh();
    utility::create_directory(output_path);

    auto content = as_pdb();
    if (content.size() == 1) {
        std::ofstream output(output_path);
        if (!output.is_open()) {throw std::ios_base::failure("PDBWriter::write: Could not open file \"" + output_path + "\"");}
        output << content[0] << std::flush;
        output.close();
        std::cout << "Output written to file " + output_path + "." << std::endl;
    }
    else {
        for (unsigned int i = 0; i < content.size(); i++) {
            std::string path = utility::stem_append(output_path, "_" + std::to_string(i+1));
            std::ofstream output(path);
            if (!output.is_open()) {throw std::ios_base::failure("PDBWriter::write: Could not open file \"" + path + "\"");}
            output << content[i] << std::flush;
            output.close();
        }
    }
}

std::vector<std::string> PDBWriter::as_pdb() const {
    File& f = *file;
    std::vector<std::string> files;
    std::string s = f.header.get();

    unsigned int count = 0;
    unsigned int i_ter = f.terminate.serial;
    bool printed_ter = false;
    for (unsigned int i = 0; i < f.protein_atoms.size(); i++) {
        if (i == i_ter) { // check if this is where the terminate is supposed to go
            s += f.terminate.as_pdb(); // write it if so
            printed_ter = true;
        }
        s += f.protein_atoms[i].as_pdb();
        count++;
        if (count == 100000) {
            count = 0;
            files.push_back(std::move(s));
            s = "";
        }
    }

    // print terminate if missing
    if (!printed_ter) {s += f.terminate.as_pdb();}

    // print hetatoms
    for (unsigned int i = 0; i < f.hydration_atoms.size(); i++) {
        s += f.hydration_atoms[i].as_pdb();
        count++;
        if (count == 100000) {
            count = 0;
            files.push_back(std::move(s));
            s = "";
        }
    }

    s += f.footer.get();
    files.push_back(std::move(s));
    return files;
}