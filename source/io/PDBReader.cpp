#include <io/PDBReader.h>
#include <io/Reader.h>
#include <io/File.h>
#include <data/Terminate.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <utility/Exceptions.h>
#include <utility/Settings.h>

#include <fstream>

void PDBReader::read(std::string input_path) {
    std::cout << "PDBReader::read" << std::endl;
    if (setting::general::verbose) {
        utility::print_info("Reading PDB file from \"" + input_path + "\"");
    }

    // check if file was succesfully opened
    std::ifstream input(input_path);
    if (!input.is_open()) {throw except::io_error("PDBReader::read: Could not open file \"" + input_path + "\"");}
    std::cout << "\tifstream successfully opened" << std::endl;

    string line; // placeholder for the current line
    File& f = *file;
    while(getline(input, line)) {
        std::cout << line << std::endl;
        string type = line.substr(0, std::min(6, int(line.size()))); // read the first 6 characters
        switch(Record::get_type(type)) {
            case Record::RecordType::ATOM: {
                std::cout << "\tRecordType::ATOM" << std::endl;
                // first just parse it as an atom; we can reuse it anyway even if it is a water molecule
                Atom atom;
                atom.parse_pdb(line);

                // if this is a water molecule, add it to the hydration atoms
                // otherwise add it to the protein atoms
                if (atom.is_water()) {
                    f.add(Water(std::move(atom)));
                } else {
                    f.add(atom);
                }
                break;
            } case Record::RecordType::TERMINATE: {
                std::cout << "\tRecordType::TERMINATE" << std::endl;
                Terminate term;
                term.parse_pdb(line);
                f.add(term);
                break;
            } case Record::RecordType::HEADER: {
                std::cout << "\tRecordType::HEADER" << std::endl;
                f.add(Record::RecordType::HEADER, line);
                break;
            } case Record::RecordType::FOOTER: {
                std::cout << "\tRecordType::FOOTER" << std::endl;
                f.add(Record::RecordType::FOOTER, line);
                break;
            } case Record::RecordType::NOTYPE: {
                std::cout << "\tRecordType::NOTYPE" << std::endl;
                break;
            } default: {
                std::cout << "\tRecordType::default" << std::endl;
                throw except::io_error("PDBReader::read: Malformed input file - unrecognized type \"" + type + "\".");
            }
        };
    }
    input.close();
    std::cout << "iterated through entire file" << std::endl;

    unsigned int n_pa = f.protein_atoms.size();
    unsigned int n_ha = f.hydration_atoms.size();
    
    if (setting::general::verbose) {
        std::cout << "\tSuccessfully read " << n_pa + n_ha << " atomic records." << std::endl;
        if (n_ha != 0) {
            std::cout << "\t\t" << f.hydration_atoms.size() << " of these are hydration atoms." << std::endl;
        }
    }
    std::cout << "PDBReader::read end" << std::endl;
}