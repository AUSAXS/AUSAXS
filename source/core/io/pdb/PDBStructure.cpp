/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/pdb/PDBStructure.h>
#include <io/pdb/Record.h>
#include <io/pdb/PDBAtom.h>
#include <io/pdb/PDBWater.h>
#include <io/ExistingFile.h>
#include <utility/Exceptions.h>
#include <data/Body.h>
#include <data/Molecule.h>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

PDBStructure::PDBStructure() = default;

PDBStructure::~PDBStructure() = default;

PDBStructure::PDBStructure(const std::vector<PDBAtom>& atoms, const std::vector<PDBWater>& waters) : atoms(atoms), waters(waters) {}

PDBStructure::PDBStructure(const std::vector<PDBAtom>& atoms, const std::vector<PDBWater>& waters, const Header& header, const Footer& footer, const Terminate& terminate) 
    : header(header), footer(footer), terminate(terminate), atoms(atoms), waters(waters) {}

auto add_single_body = [] (std::vector<PDBAtom>& atoms, std::vector<PDBWater>& waters, const data::Body& body, int& serial, int& residue_serial, char chain) {
    for (const auto& a : body.get_atoms()) {
        atoms.emplace_back(
            ++serial, form_factor::to_string(a.form_factor_type()), "", "UNK", chain, 0, "", a.coordinates(), 1, 1, form_factor::get_type(a.form_factor_type()), ""
        );
    }
    for (const auto& w : body.get_waters()) {
        waters.emplace_back(
            ++serial, "O", "", "HOH", chain, ++residue_serial, "", w.coordinates(), 1, 1, constants::atom_t::O, ""
        );
    }
};

PDBStructure::PDBStructure(const data::Body& body) {
    int serial = 0;
    int residue_serial = 0;
    atoms.reserve(body.size_atom());
    waters.reserve(body.size_water());
    add_single_body(this->atoms, this->waters, body, serial, residue_serial, 'A');
}

PDBStructure::PDBStructure(const data::Molecule& molecule) {
    int serial = 0;
    int residue_serial = 0;
    char chain = 'A';
    atoms.reserve(molecule.size_atom());
    waters.reserve(molecule.size_water());
    for (const auto& body : molecule.get_bodies()) {
        add_single_body(this->atoms, this->waters, body, serial, residue_serial, chain++);
    }
}

void PDBStructure::update(std::vector<PDBAtom>& patoms, std::vector<PDBWater>& hatoms) {
    atoms = patoms;
    waters = hatoms;
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

PDBStructure::_res PDBStructure::reduced_representation() {
    PDBStructure::_res res;
    res.atoms.reserve(atoms.size());
    res.waters.reserve(waters.size());

    for (auto& a : atoms) {
        res.atoms.emplace_back(a.coords, a.effective_charge*a.occupancy, form_factor::get_type(a.get_element(), a.get_atomic_group()));
    }

    for (auto& w : waters) {
        res.waters.emplace_back(w.coords);
    }
    return res;
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