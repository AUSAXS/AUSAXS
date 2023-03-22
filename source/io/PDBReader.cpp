#include <io/PDBReader.h>
#include <io/Reader.h>
#include <io/File.h>
#include <data/Terminate.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <utility/Exceptions.h>
#include <utility/Settings.h>
#include <utility/Constants.h>

#include <fstream>

void PDBReader::read(std::string input_path) {
    if (setting::general::verbose) {
        utility::print_info("\nReading PDB file from \"" + input_path + "\"");
    }

    // check if file was succesfully opened
    std::ifstream input(input_path);
    if (!input.is_open()) {throw except::io_error("PDBReader::read: Could not open file \"" + input_path + "\"");}

    std::string line; // placeholder for the current line
    File& f = *file;
    while(getline(input, line)) {
        std::string type = line.substr(0, std::min(6, int(line.size()))); // read the first 6 characters
        switch(Record::get_type(type)) {
            case Record::RecordType::ATOM: {
                // first just parse it as an atom; we can reuse it anyway even if it is a water molecule
                Atom atom;
                atom.parse_pdb(line);

                // if this is a water molecule, add it to the hydration atoms
                // otherwise add it to the protein atoms
                if (atom.element == constants::symbols::hydrogen) {
                    if (!setting::general::keep_hydrogens) {continue;}
                    f.add(atom);
                }
                if (atom.is_water()) {f.add(Water(std::move(atom)));} 
                else {f.add(atom);}
                break;
            } case Record::RecordType::TERMINATE: {
                Terminate term;
                term.parse_pdb(line);
                f.add(term);
                break;
            } case Record::RecordType::HEADER: {
                f.add(Record::RecordType::HEADER, line);
                break;
            } case Record::RecordType::FOOTER: {
                f.add(Record::RecordType::FOOTER, line);
                break;
            } case Record::RecordType::NOTYPE: {
                break;
            } default: {
                throw except::io_error("PDBReader::read: Malformed input file - unrecognized type \"" + type + "\".");
            }
        };
    }
    input.close();

    unsigned int n_pa = f.protein_atoms.size();
    unsigned int n_ha = f.hydration_atoms.size();
    
    if (setting::general::verbose) {
        std::cout << "\tSuccessfully read " << n_pa + n_ha << " atomic records." << std::endl;
        if (n_ha != 0) {
            std::cout << "\t\t" << f.hydration_atoms.size() << " of these are hydration atoms." << std::endl;
        }
    }
}