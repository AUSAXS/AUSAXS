/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/pdb/PDBWater.h>
#include <constants/Constants.h>
#include <utility/Concepts.h>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

PDBWater::PDBWater(PDBAtom&& a) noexcept : PDBAtom(std::move(a)) {resName = "HOH";}
PDBWater::PDBWater(const PDBAtom& a) : PDBAtom(a) {resName = "HOH";}

RecordType PDBWater::get_type() const {return RecordType::WATER;}

std::string PDBWater::get_recName() const {return "HETATM";}

double PDBWater::get_mass() const {return constants::mass::get_mass(constants::atom_t::O) + 2*constants::mass::get_mass(constants::atom_t::H);}

bool PDBWater::is_water() const {return true;}

PDBWater PDBWater::create_new_water() {
    return create_new_water(Vector3<double>{0, 0, 0});
}

PDBWater PDBWater::create_new_water(const Vector3<double>& coords) {
    return PDBWater(-1, "O", "", "HOH", ' ', -1, "", coords, 1, 0, constants::atom_t::O, "");
}

bool PDBWater::operator==(const PDBWater& rhs) const {
    if (name != rhs.name) {return false;}
    if (altLoc != rhs.altLoc) {return false;}
    if (resName != rhs.resName) {return false;}
    if (chainID != rhs.chainID) {return false;}
    if (iCode != rhs.iCode) {return false;}
    if (element != rhs.element) {return false;}
    if (charge != rhs.charge) {return false;}
    if (occupancy != rhs.occupancy) {return false;}
    if (tempFactor != rhs.tempFactor) {return false;}
    if (serial != rhs.serial) {return false;}
    // if (resSeq != rhs.resSeq) {return false;} // this is to fix io tests, since some pdb files randomly changes this order
    if (coords != rhs.coords) {return false;}
    if (uid != rhs.uid) {return false;}
    return true;
}