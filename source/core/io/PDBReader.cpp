/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/PDBReader.h>
#include <data/detail/AtomCollection.h>
#include <io/ExistingFile.h>
#include <data/record/Terminate.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <settings/GeneralSettings.h>
#include <settings/MoleculeSettings.h>
#include <constants/ConstantsFwd.h>
#include <utility/StringUtils.h>
#include <utility/Console.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::io::detail;
using namespace ausaxs::data::record;

PDBReader::PDBReader(observer_ptr<data::detail::AtomCollection> const file) : file(file) {}

PDBReader::~PDBReader() = default;

auto parse_single_file = [] (const io::ExistingFile& file, data::detail::AtomCollection& collection) -> void {
    // check if file was succesfully opened
    std::ifstream input(file);
    if (!input.is_open()) {throw except::io_error("PDBReader::read: Could not open file \"" + file + "\"");}

    unsigned int discarded_hydrogens = 0;
    std::string line;
    while(getline(input, line)) {
        if (utility::remove_all(line, " \n\r").empty()) {continue;}
        std::string type = line.substr(0, std::min(6, int(line.size()))); // read the first 6 characters
        switch(Record::get_type(type)) {
            case RecordType::ATOM: {
                // first just parse it as an atom; we can reuse it anyway even if it is a water molecule
                Atom atom;
                atom.parse_pdb(line);

                // check if this is a hydrogen atom
                if (atom.element == constants::atom_t::H && !settings::general::keep_hydrogens) {
                    discarded_hydrogens++;
                    continue;
                }

                // check if this is a water molecule
                if (atom.is_water()) {collection.add(Water(std::move(atom)));} 
                else {collection.add(std::move(atom));}
                break;
            } case RecordType::TERMINATE: {
                Terminate term;
                term.parse_pdb(line);
                collection.add(term);
                break;
            } case RecordType::HEADER: {
                collection.add(RecordType::HEADER, line);
                break;
            } case RecordType::FOOTER: {
                collection.add(RecordType::FOOTER, line);
                break;
            } case RecordType::NOTYPE: {
                break;
            } default: {
                throw except::io_error("PDBReader::read: Malformed input file - unrecognized type \"" + type + "\".");
            }
        };
    }
    input.close();
    
    if (!settings::molecule::use_occupancy) {
        for (auto& a : collection.protein_atoms) {a.occupancy = 1.0;}
    }

    if (discarded_hydrogens != 0) {
        console::print_text("Discarded " + std::to_string(discarded_hydrogens) + " explicit hydrogen atoms.");
    }
};

void PDBReader::read(const io::File& path) {
    console::print_info("Reading PDB file from \"" + path.str() + "\"");
    console::indent();

    if (path.append("_part1").exists()) {
        console::print_text("File is split into multiple parts.");
        unsigned int i = 1;
        while (path.append("_part" + std::to_string(i)).exists()) {
            console::print_text("\tParsed file " + path.append("_part" + std::to_string(i)).str());
            parse_single_file(path.append("_part" + std::to_string(i)), *file);
            i++;
        }
    } else {
        parse_single_file(path, *file);
    }

    unsigned int n_pa = file->protein_atoms.size();
    unsigned int n_ha = file->hydration_atoms.size();

    console::print_text("Successfully read " + std::to_string(n_pa + n_ha) + " atomic records.");
    if (n_ha != 0) {console::print_text("\t" + std::to_string(file->hydration_atoms.size()) + " of these are hydration atoms.");}
    console::unindent();
}