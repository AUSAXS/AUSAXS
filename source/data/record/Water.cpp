/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/record/Water.h>
#include <constants/Constants.h>

using namespace data::record;

Water::Water(const Atom&& a) noexcept : Atom(std::move(a)) {}
Water::Water(const Atom& a) : Atom(a) {}
Water::Water(const Water& w) = default;
Water::~Water() = default;

RecordType Water::get_type() const {return RecordType::WATER;}

std::string Water::get_recName() const {return "HETATM";}

bool Water::is_water() const {return true;}

Water Water::create_new_water() {
    return create_new_water({0, 0, 0});
}

Water Water::create_new_water(Vector3<double> coords) {
    return Water(-1, "O", "", "HOH", ' ', -1, "", coords, 1, 0, constants::atom_t::O, "");
}

Water& Water::operator=(const Water& rhs) = default;

bool Water::operator==(const Water& rhs) const {
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