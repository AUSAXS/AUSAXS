/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/Body.h>
#include <data/BodySymmetryFacade.h>
#include <data/state/UnboundSignaller.h>
#include <data/Symmetry.h>
#include <grid/Grid.h>
#include <constants/Constants.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <math/Vector3.h>

#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::data;

Body::Body() : uid(uid_counter++) {initialize();}
Body::Body(const Body& body) : uid(body.uid) {initialize();}
Body::Body(Body&& body) noexcept : uid(body.uid) {initialize();}
Body::~Body() = default;

Body::Body(const io::File& path) : uid(uid_counter++) {
    initialize();
}

Body::Body(const std::vector<data::AtomFF>& atoms, const std::vector<data::Water>& waters) : uid(uid_counter++) {
    initialize();
}

Body::Body(const std::vector<data::AtomFF>& atoms) : Body(atoms, std::vector<data::Water>()) {}

void Body::initialize() {
    signal = std::make_shared<signaller::UnboundSignaller>();
}

data::detail::BodySymmetryFacade<Body> Body::symmetry() {return data::detail::BodySymmetryFacade<Body>(this);}

data::detail::BodySymmetryFacade<const Body> Body::symmetry() const {return data::detail::BodySymmetryFacade<const Body>(this);}

Vector3<double> Body::get_cm() const {
    Vector3<double> cm{0, 0, 0};
    double M = 0; // total mass
    auto weighted_sum = [&cm, &M] (auto& atoms) {
        for (auto const& a : atoms) {
            double m = constants::mass::get_mass(a.form_factor_type());
            M += m;
            cm += a.coordinates()*m;
        }
    };
    weighted_sum(atoms);
    weighted_sum(waters);
    return cm/M;
}

double Body::get_volume_vdw() const {
    double volume = std::accumulate(atoms.begin(), atoms.end(), 0.0, [] (double sum, const data::AtomFF& atom) {
        return sum + std::pow(constants::radius::get_vdw_radius(atom.form_factor_type()), 3);
    });
    return 4*std::numbers::pi*volume/3;
}

void Body::translate(Vector3<double> v) {
    signal->external_change();
    std::for_each(atoms.begin(), atoms.end(), [v] (data::AtomFF& atom) {atom.coordinates() += v;});
    std::for_each(waters.begin(), waters.end(), [v] (data::Water& atom) {atom.coordinates() += v;});
}

void Body::rotate(const Matrix<double>& R) {
    signal->external_change();
    for (auto& atom : atoms) {
        atom.coordinates().rotate(R);
    }

    for (auto& atom : waters) {
        atom.coordinates().rotate(R);
    }
}

double Body::get_total_atomic_charge() const {
    return std::accumulate(get_atoms().begin(), get_atoms().end(), 0.0, [] (double sum, const data::AtomFF& atom) {return sum + atom.weight();});
}

double Body::get_molar_mass() const {
    return get_absolute_mass()*constants::Avogadro;
}

double Body::get_absolute_mass() const {
    double M = 0;
    std::for_each(atoms.begin(), atoms.end(), [&M] (const data::AtomFF& a) {M += constants::mass::get_mass(a.form_factor_type());});
    std::for_each(waters.begin(), waters.end(), [&M] (const data::Water& a) {M += constants::mass::get_mass(a.form_factor_type());});
    return M;
}

Body& Body::operator=(Body&& rhs) noexcept {
    atoms = std::move(rhs.atoms);
    waters = std::move(rhs.waters);
    symmetries = std::move(rhs.symmetries);
    uid = rhs.uid;    
    signal->internal_change();
    signal->external_change();
    return *this;
}

Body& Body::operator=(const Body& rhs) {
    atoms = std::move(rhs.atoms);
    waters = std::move(rhs.waters);
    symmetries = std::move(rhs.symmetries);
    signal->internal_change();
    signal->external_change();
    return *this;
}

bool Body::operator==(const Body& rhs) const {
    return uid == rhs.uid;
}

bool Body::equals_content(const Body& rhs) const {
    return atoms == rhs.atoms && waters == rhs.waters && symmetries == rhs.symmetries;
}

std::shared_ptr<signaller::Signaller> Body::get_signaller() const {
    return signal;
}

void Body::register_probe(std::shared_ptr<signaller::Signaller> signal) {
    this->signal = std::move(signal);
}

std::vector<data::AtomFF>& Body::get_atoms() {return atoms;}

std::vector<data::Water>& Body::get_waters() {return waters;}

const std::vector<data::AtomFF>& Body::get_atoms() const {return atoms;}

const std::vector<data::Water>& Body::get_waters() const {return waters;}

data::AtomFF& Body::get_atom(unsigned int index) {return atoms[index];}

const data::AtomFF& Body::get_atom(unsigned int index) const {return atoms[index];}

int Body::get_uid() const {return uid;}

std::size_t Body::size_atom() const {return atoms.size();}

std::size_t Body::size_water() const {return waters.size();}

std::size_t Body::size_symmetry() const {return symmetries.size();}

std::size_t Body::size_symmetry_total() const {
    return std::accumulate(symmetries.begin(), symmetries.end(), 0, [] (int sum, const detail::Symmetry& sym) {return sum + sym.repeat;});
}