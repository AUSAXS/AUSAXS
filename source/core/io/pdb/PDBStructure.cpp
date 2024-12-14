/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/pdb/PDBStructure.h>
#include <io/pdb/Record.h>
#include <io/pdb/PDBAtom.h>
#include <io/pdb/PDBWater.h>
#include <io/pdb/Terminate.h>
#include <io/pdb/Header.h>
#include <io/pdb/Footer.h>
#include <io/ExistingFile.h>
#include <io/PDBReader.h>
#include <io/CIFReader.h>
#include <io/PDBWriter.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

PDBStructure::PDBStructure() = default;

PDBStructure::PDBStructure(const std::vector<PDBAtom>& atoms, const std::vector<PDBWater>& waters) : atoms(atoms), waters(waters) {}

PDBStructure::PDBStructure(const std::vector<PDBAtom>& atoms, const std::vector<PDBWater>& waters, const Header& header, const Footer& footer, const Terminate& terminate) 
    : header(header), footer(footer), terminate(terminate), atoms(atoms), waters(waters) {}

PDBStructure::PDBStructure(const io::File& file) {
    read(file);
}

void PDBStructure::update(std::vector<PDBAtom>& patoms, std::vector<PDBWater>& hatoms) {
    atoms = patoms;
    waters = hatoms;
}

void PDBStructure::read(const io::File& path) {
    construct_reader(path)->read(path);
}

void PDBStructure::write(const io::File& path) {
    refresh();
    construct_writer(path)->write(path);
}

const std::vector<PDBAtom>& PDBStructure::get_atoms() const {return atoms;}

const std::vector<PDBWater> PDBStructure::get_waters() const {return waters;}

void PDBStructure::add(const PDBAtom& r) {
    atoms.push_back(r);
}

void PDBStructure::add(PDBAtom&& r) {
    atoms.push_back(std::move(r));
}

void PDBStructure::add(const PDBWater& r) {
    waters.push_back(r);
}

void PDBStructure::add(PDBWater&& r) {
    waters.push_back(std::move(r));
}

void PDBStructure::add(const Terminate& ter) {
    terminate = ter;
    // terminates.push_back(r);
}

void PDBStructure::add(const RecordType& type, const std::string& s) {
    if (type == RecordType::HEADER) {
        header.add(s);
    } else if (type == RecordType::FOOTER) {
        footer.add(s);
    } else {
        throw except::invalid_argument("PDBStructure::add: Type is not \"HEADER\" or \"FOOTER\"!");
    }
}

void PDBStructure::refresh() {
    if (atoms.empty()) {
        terminate = Terminate(-1, "", ' ', -1, ""); 
        return;
    }

    bool terminate_inserted = false;
    char chainID = '0'; int resSeq = 0; int serial = atoms[0].serial;

    auto insert_ter = [&] () {
        // last atom before the terminate
        // we need this to determine what chainID and resSeq to use for the terminate and hetatms
        const PDBAtom& a = atoms.at(serial-1-atoms[0].serial);
        chainID = a.chainID;
        resSeq = a.resSeq;
        if (serial != 0) {terminate = Terminate(serial, a.resName, a.chainID, a.resSeq, " ");}
        terminate_inserted = true;
    };

    for (auto& a : atoms) {
        if (!terminate_inserted && a.get_type() == RecordType::WATER) {
            insert_ter();
            resSeq++; // TER records always denotes the end of a sequence
            serial++;
        }
        a.set_serial(serial++ % 100000); // fix possible errors in the serial
    }

    if (!terminate_inserted) {
        insert_ter();
        resSeq++; // TER records always denotes the end of a sequence
        serial++;
    }

    chainID = atoms[atoms.size()-1].chainID+1;
    resSeq = atoms[atoms.size()-1].resSeq + 1;
    for (auto& a : waters) {
        a.set_serial(serial++ % 100000);
        a.set_residue_sequence_number(resSeq++ % 10000);
        a.set_chainID(chainID + int(resSeq/10000));
    }
}

bool PDBStructure::operator==(const PDBStructure& rhs) const = default;

#define FAILURE_MSG false
#if FAILURE_MSG
    #include <iostream>
#endif
bool PDBStructure::equals_content(const PDBStructure& rhs) const {
    if (atoms.size() != rhs.atoms.size()) {
        #if FAILURE_MSG
            std::cout << "atoms.size() != rhs.atoms.size()" << std::endl;
        #endif
        return false;
    }

    if (waters.size() != rhs.waters.size()) {
        #if FAILURE_MSG
            std::cout << "waters.size() != rhs.waters.size()" << std::endl;
        #endif
        return false;
    }

    for (unsigned int i = 0; i < atoms.size(); i++) {
        if (!atoms[i].equals_content(rhs.atoms[i])) {
            #if FAILURE_MSG
                std::cout << "!atoms[" << i << "].equals_content(rhs.atoms[" << i << "])" << std::endl;
            #endif
            return false;
        }
    }

    for (unsigned int i = 0; i < waters.size(); i++) {
        if (!waters[i].equals_content(rhs.waters[i])) {
            #if FAILURE_MSG
                std::cout << "!waters[" << i << "].equals_content(rhs.waters[" << i << "])" << std::endl;
            #endif
            return false;
        }
    }

    return true;    
}

std::unique_ptr<io::detail::Reader> PDBStructure::construct_reader(const io::File& path) {
    if (path.extension() == ".xml" || path.extension() == ".XML") { // .xml PDBStructure
        throw except::invalid_argument("PDBStructure::construct_reader: .xml input PDBStructures are not supported.");
    } else if (path.extension() == ".pdb" || path.extension() == ".PDB") { // .pdb PDBStructure
        return std::make_unique<io::detail::PDBReader>(this);
    } else if (path.extension() == ".ent" || path.extension() == ".ENT") { // .ent PDBStructure
        return std::make_unique<io::detail::PDBReader>(this);
    } else if (path.extension() == ".cif" || path.extension() == ".CIF") { // .cif PDBStructure
        return std::make_unique<io::detail::CIFReader>(this);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("PDBStructure::construct_reader: Unsupported extension \"" + path.extension() + "\" of input file \"" + path.str() + "\".");
    }
}

std::unique_ptr<io::detail::Writer> PDBStructure::construct_writer(const io::File& path) {
    if (path.extension() == ".xml" || path.extension() == ".XML") { // .xml PDBStructure
        throw except::invalid_argument("PDBStructure::construct_writer: .xml input PDBStructures are not supported.");
    } else if (path.extension() == ".pdb" || path.extension() == ".PDB") { // .pdb PDBStructure
        return std::make_unique<io::detail::PDBWriter>(this);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("PDBStructure::construct_writer: Unsupported extension \"" + path.extension() + "\" of input file \"" + path.str() + "\".");
    }
}