// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#ifdef _MSC_VER
    #pragma warning(disable:4996) // disable sscanf deprecation warning on MSVC
#endif

#include <io/pdb/PDBAtom.h>
#include <constants/Constants.h>
#include <utility/Utility.h>
#include <settings/MoleculeSettings.h>
#include <utility/Console.h>

#include <utility>
#include <iomanip>
#include <iostream>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

PDBAtom::PDBAtom() : uid(uid_counter++) {}

PDBAtom::PDBAtom(Vector3<double> v, double occupancy, constants::atom_t element, const std::string& resName, int serial) : uid(uid_counter++) {
    // we use our setters so we can validate the input if necessary
    this->coords = v;
    this->occupancy = occupancy;
    this->element = element;
    this->resName = resName;
    this->serial = serial;
    this->effective_charge = constants::charge::nuclear::get_charge(this->element);
}

PDBAtom::PDBAtom(int serial, const std::string& name, const std::string& altLoc, const std::string& resName, char chainID, int resSeq, const std::string& iCode, 
    Vector3<double> coords, double occupancy, double tempFactor, constants::atom_t element, const std::string& charge) : uid(uid_counter++) 
{
    this->serial = serial;
    this->name = name;
    this->altLoc = altLoc;
    this->resName = resName;
    this->chainID = chainID;
    this->resSeq = resSeq;
    this->iCode = iCode;
    this->coords = coords;
    this->occupancy = occupancy;
    this->tempFactor = tempFactor;
    this->element = element;
    this->charge = charge;
    this->effective_charge = constants::charge::nuclear::get_charge(this->element);
    atomic_group = constants::atomic_group_t::unknown;
}

void PDBAtom::parse_pdb(const std::string& str) {
    auto s = utility::remove_all(str, "\n\r"); // remove any newline or carriage return
    int pad_size = 80 - static_cast<int>(s.size());
    if (pad_size < 0) {
        static bool warned = false;
        if (!warned) {
            console::print_warning("Warning in PDBAtom::parse_pdb: Found line longer than 80 characters. Truncating. Further warnings of this type will be suppressed.");
            warned = true;
        }
        s = s.substr(0, 80);
    } else if (pad_size > 0) {
        static bool warned = false;
        if (!warned) {
            console::print_warning("Warning in PDBAtom::parse_pdb: Found line shorter than 80 characters. Padding with spaces. Further warnings of this type will be suppressed.");
            warned = true;
        }
        s += std::string(pad_size, ' ');
    }

    // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    //                   RN SE S1 NA AL RN CI S2 RS iC S3 X  Y  Z  OC TF S4  EL CH
    //                   0     1           2              3     4  5  6      7     8
    //                   0  6  1  2  6  7  0  1  2  6  7  0  8  6  4  0  6   6  8  0  
    std::string form = "%6c%5c%1c%4c%1c%3c%1c%1c%4c%1c%3c%8c%8c%8c%6c%6c%10c%2c%2c";
    std::string recName = "      ", serial = "     ", space1 = " ", name = "    ", altLoc = " ", resName = "   ", space2 = " ", 
        resSeq = "    ", iCode = " ", space3 = "   ", x = "        ", y = "        ", z = "        ", 
        occupancy = "      ", tempFactor = "      ", space4 = "          ", element = "  ", charge = "  ";
    char chainID = ' ';
    sscanf(s.c_str(), form.data(), recName.data(), serial.data(), space1.data(), name.data(), altLoc.data(), resName.data(), 
        space2.data(), &chainID, resSeq.data(), iCode.data(), space3.data(), x.data(), y.data(), z.data(), 
        occupancy.data(), tempFactor.data(), space4.data(), element.data(), charge.data());

    // sanity check
    if (!(Record::get_type(recName) == RecordType::ATOM)) [[unlikely]] {
        throw except::parse_error("PDBAtom::parse_pdb: input std::string is not \"ATOM  \" or \"HETATM\" (" + recName + ").");
    }

    // remove any spaces from the numbers
    serial = utility::remove_all(serial, " ");
    name = utility::remove_all(space1+name, " "); // we add space1 since some programs (gromacs) uses it for the name.
    resName = utility::remove_all(resName, " ");
    resSeq = utility::remove_all(resSeq, " ");
    x = utility::remove_all(x, " ");
    y = utility::remove_all(y, " ");
    z = utility::remove_all(z, " ");
    occupancy = utility::remove_all(occupancy, " ");
    tempFactor = utility::remove_all(tempFactor, " ");
    element = utility::remove_all(element, " ");

    // sometimes people use the first character of x for some other purpose.
    // if it is a digit, the following won't work. On the other hand they're kinda asking for it then. Follow the standard, people. 
    if (!(std::isdigit(x[0]) || x[0] == '-')) {
        x = x.substr(1, x.size()-1);
    }

    // set all of the properties
    try {
        this->serial = std::stoi(serial);
        this->name = std::move(name);
        this->altLoc = std::move(altLoc);
        this->resName = std::move(resName);
        this->chainID = chainID;
        this->resSeq = std::stoi(resSeq);
        this->iCode = std::move(iCode);
        this->coords = {std::stod(x), std::stod(y), std::stod(z)};
        if (occupancy.empty()) {this->occupancy = 1;} else {this->occupancy = std::stod(occupancy);}
        if (tempFactor.empty()) {this->tempFactor = 0;} else {this->tempFactor = std::stod(tempFactor);}
        if (element.empty()) {
            // if the element is not set, we can try to infer it from the name
            if (auto s = this->name.substr(0, 1); !std::isdigit(s[0])) [[likely]] {
                set_element(s);
            } else {
                set_element(this->name.substr(1, 1)); // sometimes the first character is a number
            }
        } else {
            set_element(element);
        }
        this->charge = std::move(charge);
    } catch (const except::base& e) { // catch conversion errors and output a more meaningful error message
        console::print_warning("PDBAtom::parse_pdb: Invalid field values in line \"" + s + "\".");
        throw e;
    }

    effective_charge = constants::charge::nuclear::get_charge(this->element);
    atomic_group = constants::atomic_group_t::unknown;
}

void PDBAtom::add_implicit_hydrogens() {
    assert(element != constants::atom_t::H && "PDBAtom::add_implicit_hydrogens: Attempted to add implicit hydrogens to a hydrogen atom.");
    try {
        effective_charge = constants::charge::nuclear::get_charge(element) + constants::hydrogen_atoms::residues.get(resName).get(name, element);
        atomic_group = constants::symbols::get_atomic_group(resName, name, element);
    } catch (const except::base&) {
        throw except::invalid_argument(
            "PDBAtom::add_implicit_hydrogens: Could not identify group of atom " + std::to_string(serial) + ". Unknown element, residual or atom: "
            "(" + constants::symbols::to_string(element) + ", " + this->resName + ", " + this->name + ")"
        );
    }
}

using std::left, std::right, std::setw;
std::string PDBAtom::as_pdb() const {
    std::stringstream ss;
    //                   RN SE S1 NA AL RN CI S2 RS iC S3 X  Y  Z  OC TF S2  EL CH
    //                   0     1           2              3     4  5  6      7     8
    //                   0  6  1  2  6  7  0  1  2  6  7  0  8  6  4  0  6   6  8  0  
    //          format: "%6c%5c%2c%4c%1c%3c %1c%4c%1c%3c%8c%8c%8c%6c%6c%10c%2c%2c"
    ss << left << setw(6) << get_recName()                                          // 1 - 6
        << right << setw(5) << serial                                               // 7 - 11
        << " "                                                                      // 12
        << " " << left << setw(3) << name                                           // 13 - 16
        << left << setw(1) << altLoc                                                // 17
        << left << setw(3) << resName                                               // 18 - 20
        << " "                                                                      // 21
        << left << setw(1) << chainID                                               // 22
        << right << setw(4) << resSeq                                               // 23 - 26
        << right << setw(1) << iCode                                                // 27
        << "   "                                                                    // 28 - 30
        << right << setw(8) << utility::fixedwidth(coords.x(), 7)                   // 31 - 38
        << right << setw(8) << utility::fixedwidth(coords.y(), 7)                   // 39 - 46
        << right << setw(8) << utility::fixedwidth(coords.z(), 7)                   // 47 - 54
        << right << setw(6) << utility::fixedwidth(occupancy, 6)                    // 55 - 60
        << right << setw(6) << utility::fixedwidth(tempFactor, 6)                   // 61 - 66
        << "          "                                                             // 67 - 76
        << right << setw(2) << constants::symbols::to_string(element)               // 77 - 78
        << left << setw(2) << charge                                                // 79 - 80
        << std::endl;
    return ss.str();
}

std::string PDBAtom::get_recName() const {return "ATOM  ";}

RecordType PDBAtom::get_type() const {return RecordType::ATOM;}

bool PDBAtom::is_water() const {return (resName == "HOH") || (resName == "SOL");}

void PDBAtom::set_element(constants::atom_t element) {
    assert(element != constants::atom_t::unknown && "PDBAtom::set_element: Attempted to set element to unknown.");
    this->element = element;
}

void PDBAtom::set_element(const std::string& element) {
    set_element(constants::symbols::parse_element_string(element));
}

Vector3<double>& PDBAtom::coordinates() {return coords;}
const Vector3<double>& PDBAtom::coordinates() const {return coords;}

double PDBAtom::get_mass() const {
    if (settings::molecule::implicit_hydrogens) {
        #ifdef DEBUG
            try {
                return constants::mass::get_mass(element) + constants::hydrogen_atoms::residues.get(this->resName).get(this->name, this->element)*constants::mass::get_mass(constants::atom_t::H);
            } catch (const std::exception& e) {
                console::print_critical(e.what());
                throw except::invalid_argument("PDBAtom::get_mass: The mass of element " + constants::symbols::to_string(element) + " (serial " + std::to_string(serial) + ") is not defined.");
            }
        #endif
        // mass of this nucleus + mass of attached H atoms
        return constants::mass::get_mass(element) + constants::hydrogen_atoms::residues.get(this->resName).get(this->name, this->element)*constants::mass::get_mass(constants::atom_t::H);
    } else {
        #ifdef DEBUG
            if (element == constants::atom_t::unknown) [[unlikely]] {
                throw except::invalid_argument("PDBAtom::get_mass: Attempted to get atomic mass, but the element was not set!");
            }
        #endif
        return constants::mass::get_mass(element);
    }
}

unsigned int PDBAtom::Z() const {
    #ifdef DEBUG
        if (element == constants::atom_t::unknown) [[unlikely]] {
            throw except::invalid_argument("PDBAtom::get_Z: Attempted to get atomic charge, but the element was not set!");
        }
    #endif
    return constants::charge::nuclear::get_charge(element);
}

bool PDBAtom::operator<(const PDBAtom& rhs) const {
    return serial < rhs.serial;
}

bool PDBAtom::operator==(const PDBAtom& rhs) const {
    return uid == rhs.uid;
}

#define FAILURE_MSG false
#if FAILURE_MSG
    #include <iostream>
#endif
bool PDBAtom::equals_content(const PDBAtom& rhs) const {
    if (coords != rhs.coords) {
        #if FAILURE_MSG
            std::cout << "coords \"" << coords.x() << ", " << coords.y() << ", " << coords.z() <<  "\" != rhs.coords \"" << rhs.coords.x() << ", " << rhs.coords.y() << ", " << rhs.coords.z() << "\"" << std::endl;
        #endif
        return false;
    }

    if (name != rhs.name) {
        #if FAILURE_MSG
            std::cout << "name \"" + name +  "\" != rhs.name \"" + rhs.name + "\"" << std::endl;
        #endif
        return false;
    }

    if (altLoc != rhs.altLoc) {
        #if FAILURE_MSG
            std::cout << "altLoc \"" + altLoc +  "\" != rhs.altLoc \"" + rhs.altLoc + "\"" << std::endl;
        #endif
        return false;
    }

    if (resName != rhs.resName) {
        #if FAILURE_MSG
            std::cout << "resName \"" + resName +  "\" != rhs.resName \"" + rhs.resName + "\"" << std::endl;
        #endif
        return false;
    }

    if (chainID != rhs.chainID) {
        #if FAILURE_MSG
            std::cout << "chainID \"" + std::string(1, chainID) +  "\" != rhs.chainID \"" + std::string(1, rhs.chainID) + "\"" << std::endl;
        #endif
        return false;
    }

    if (iCode != rhs.iCode) {
        #if FAILURE_MSG
            std::cout << "iCode \"" + iCode +  "\" != rhs.iCode \"" + rhs.iCode + "\"" << std::endl;
        #endif
        return false;
    }

    if (element != rhs.element) {
        #if FAILURE_MSG
            std::cout << "element \"" + constants::symbols::to_string(element) +  "\" != rhs.element \"" + constants::symbols::to_string(rhs.element) + "\"" << std::endl;
        #endif
        return false;
    }

    if (charge != rhs.charge) {
        #if FAILURE_MSG
            std::cout << "charge \"" + charge +  "\" != rhs.charge \"" + rhs.charge + "\"" << std::endl;
        #endif
        return false;
    }

    if (occupancy != rhs.occupancy) {
        #if FAILURE_MSG
            std::cout << "occupancy \"" + std::to_string(occupancy) +  "\" != rhs.occupancy \"" + std::to_string(rhs.occupancy) + "\"" << std::endl;
        #endif
        return false;
    }

    if (tempFactor != rhs.tempFactor) {
        #if FAILURE_MSG
            std::cout << "tempFactor \"" + std::to_string(tempFactor) +  "\" != rhs.tempFactor \"" + std::to_string(rhs.tempFactor) + "\"" << std::endl;
        #endif
        return false;
    }

    if (serial != rhs.serial) {
        #if FAILURE_MSG
            std::cout << "serial \"" + std::to_string(serial) +  "\" != rhs.serial \"" + std::to_string(rhs.serial) + "\"" << std::endl;
        #endif
        return false;
    }

    if (resSeq != rhs.resSeq) {
        #if FAILURE_MSG
            std::cout << "resSeq \"" + std::to_string(resSeq) +  "\" != rhs.resSeq \"" + std::to_string(rhs.resSeq) + "\"" << std::endl;
        #endif
        return false;
    }

    if (effective_charge != rhs.effective_charge) {
        #if FAILURE_MSG
            std::cout << "effective_charge \"" + std::to_string(effective_charge) +  "\" != rhs.effective_charge \"" + std::to_string(rhs.effective_charge) + "\"" << std::endl;
        #endif
        return false;
    }

    return true;
}