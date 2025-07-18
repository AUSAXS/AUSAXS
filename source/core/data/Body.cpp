// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/Body.h>
#include <data/state/UnboundSignaller.h>
#include <grid/Grid.h>
#include <constants/Constants.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <math/Vector3.h>
#include <hydrate/ExplicitHydration.h>
#include <hydrate/ImplicitHydration.h>
#include <io/Reader.h>
#include <settings/MoleculeSettings.h>

#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::data;

Body::Body() : hydration(hydrate::Hydration::create()), symmetries(std::make_unique<symmetry::SymmetryStorage>()), uid(uid_counter++) {initialize();}
Body::Body(const Body& body) : atoms(body.atoms), hydration(body.hydration->clone()), symmetries(body.symmetries->clone()), uid(body.uid) {initialize();}
Body::Body(Body&& body) noexcept : atoms(std::move(body.atoms)), hydration(std::move(body.hydration)), symmetries(std::move(body.symmetries)), uid(body.uid) {initialize();}
Body::~Body() = default;

Body::Body(const io::File& path) : uid(uid_counter++) {
    auto file = io::Reader::read(path);
    if (settings::molecule::implicit_hydrogens) {file.add_implicit_hydrogens();}
    auto data = file.reduced_representation();
    atoms = std::move(data.atoms);
    if (data.waters.empty()) {
        hydration = hydrate::Hydration::create();
    } else {
        hydration = hydrate::Hydration::create(std::move(data.waters));
    }
    symmetries = std::make_unique<symmetry::SymmetryStorage>();
    initialize();
}

auto convert_atom_atomff = [] (const std::vector<data::Atom>& atoms) {
    std::vector<data::AtomFF> atoms_ff;
    atoms_ff.reserve(atoms.size());
    for (auto& a : atoms) {
        atoms_ff.emplace_back(a, form_factor::form_factor_t::UNKNOWN);
    }
    return atoms_ff;
};

template<AtomVectorFF T>
Body::Body(T&& atoms) 
    : atoms(std::forward<T>(atoms)), hydration(hydrate::Hydration::create()), 
      symmetries(std::make_unique<symmetry::SymmetryStorage>()), uid(uid_counter++
) {
    initialize();
}

template<AtomVectorFF T, WaterVector U>
Body::Body(T&& atoms, U&& waters) 
    : atoms(std::forward<T>(atoms)), hydration(hydrate::Hydration::create(std::forward<U>(waters))), 
      symmetries(std::make_unique<symmetry::SymmetryStorage>()), uid(uid_counter++
) {
    initialize();
}

Body::Body(const std::vector<Atom>& atoms) 
    : atoms(convert_atom_atomff(atoms)), hydration(hydrate::Hydration::create()), 
      symmetries(std::make_unique<symmetry::SymmetryStorage>()), uid(uid_counter++
) {
    initialize();
}

template<WaterVector U>
Body::Body(const std::vector<Atom>& atoms, U&& waters) 
    : atoms(convert_atom_atomff(atoms)), hydration(hydrate::Hydration::create(std::forward<U>(waters))), 
    symmetries(std::make_unique<symmetry::SymmetryStorage>()), uid(uid_counter++
) {
    initialize();
}

template Body::Body(std::vector<data::AtomFF>&& atoms);
template Body::Body(const std::vector<data::AtomFF>& atoms);
template Body::Body(std::vector<data::AtomFF>& atoms);

template Body::Body(std::vector<data::AtomFF>&& atoms, std::vector<data::Water>&& waters);
template Body::Body(const std::vector<data::AtomFF>& atoms, std::vector<data::Water>&& waters);
template Body::Body(std::vector<data::AtomFF>& atoms, std::vector<data::Water>&& waters);
template Body::Body(std::vector<data::AtomFF>&& atoms, const std::vector<data::Water>& waters);
template Body::Body(std::vector<data::AtomFF>&& atoms, std::vector<data::Water>& waters);
template Body::Body(std::vector<data::AtomFF>& atoms, std::vector<data::Water>& waters);
template Body::Body(const std::vector<data::AtomFF>& atoms, const std::vector<data::Water>& waters);

template Body::Body(const std::vector<data::Atom>& atoms, std::vector<data::Water>&& waters);
template Body::Body(const std::vector<data::Atom>& atoms, const std::vector<data::Water>& waters);
template Body::Body(const std::vector<data::Atom>& atoms, std::vector<data::Water>& waters);

void Body::initialize() {
    signal = std::make_shared<signaller::UnboundSignaller>();
}

symmetry::detail::BodySymmetryFacade<Body> Body::symmetry() {return symmetry::detail::BodySymmetryFacade<Body>(this);}

symmetry::detail::BodySymmetryFacade<const Body> Body::symmetry() const {return symmetry::detail::BodySymmetryFacade<const Body>(this);}

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
    if (auto h = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get()); h) {
        weighted_sum(h->waters);
    }
    return cm/M;
}

double Body::get_volume_vdw() const {
    double volume = std::accumulate(atoms.begin(), atoms.end(), 0.0, [] (double sum, const data::AtomFF& atom) {
        return sum + std::pow(constants::radius::get_vdw_radius(atom.form_factor_type()), 3);
    });
    return 4*std::numbers::pi*volume/3;
}

void Body::translate(Vector3<double> v) {
    signal->modified_external();
    std::for_each(atoms.begin(), atoms.end(), [v] (data::AtomFF& atom) {atom.coordinates() += v;});
    if (auto h = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get()); h) {
        signal->modified_hydration();
        std::for_each(h->waters.begin(), h->waters.end(), [v] (data::Water& atom) {atom.coordinates() += v;});
    }
}

void Body::rotate(const Matrix<double>& R) {
    signal->modified_external();
    for (auto& atom : atoms) {
        atom.coordinates().rotate(R);
    }

    if (auto h = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get()); h) {
        signal->modified_hydration();
        for (auto& atom : h->waters) {
            atom.coordinates().rotate(R);
        }
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
    if (auto h = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get()); h) {
        std::for_each(h->waters.begin(), h->waters.end(), [&M] (const data::Water& a) {M += constants::mass::get_mass(a.form_factor_type());});
    }
    return M;
}

Body& Body::operator=(Body&& rhs) noexcept {
    atoms = std::move(rhs.atoms);
    hydration = std::move(rhs.hydration);
    symmetries = std::move(rhs.symmetries);
    uid = rhs.uid;
    signal->modified_internal();
    signal->modified_external();
    signal->modified_hydration();
    return *this;
}

Body& Body::operator=(const Body& rhs) {
    atoms = rhs.atoms;
    if (auto h = dynamic_cast<hydrate::ExplicitHydration*>(rhs.hydration.get()); h) {
        hydration = std::make_unique<hydrate::ExplicitHydration>(h->waters);
    } else if (auto h = dynamic_cast<hydrate::ImplicitHydration*>(rhs.hydration.get()); h) {
        throw std::runtime_error("Body::operator=: Implicit hydration is not implemented.");
    }
    symmetries = rhs.symmetries->clone();
    signal->modified_internal();
    signal->modified_external();
    signal->modified_hydration();
    return *this;
}

bool Body::operator==(const Body& rhs) const {
    return uid == rhs.uid;
}

bool Body::equals_content(const Body& rhs) const {
    if (atoms != rhs.atoms) {return false;}
    if (auto h = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get()); h) {
        if (auto r = dynamic_cast<hydrate::ExplicitHydration*>(rhs.hydration.get()); r) {
            if (h->waters != r->waters) {
                return false;
            }
        } else {
            return false;
        }
    } else if (auto h = dynamic_cast<hydrate::ImplicitHydration*>(hydration.get()); h) {
        if (auto r = dynamic_cast<hydrate::ImplicitHydration*>(rhs.hydration.get()); r) {
            throw std::runtime_error("Body::equals_content: Implicit hydration is not implemented.");
        } else {
            return false;
        }
    }
    if (symmetries->get() != rhs.symmetries->get()) {return false;}
    return true;
}

std::shared_ptr<signaller::Signaller> Body::get_signaller() const {
    return signal;
}

void Body::register_probe(std::shared_ptr<signaller::Signaller> signal) {
    this->signal = std::move(signal);
}

std::vector<data::AtomFF>& Body::get_atoms() {return atoms;}

std::optional<std::reference_wrapper<const std::vector<data::Water>>> Body::get_waters() const {
    assert(hydration != nullptr && "Body::get_waters: hydration is nullptr.");
    auto h = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get());
    return h ? std::optional(std::cref(h->waters)) : std::nullopt;
}

std::optional<std::reference_wrapper<std::vector<data::Water>>> Body::get_waters() {
    assert(hydration != nullptr && "Body::get_waters: hydration is nullptr.");
    auto h = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get());
    return h ? std::optional(std::ref(h->waters)) : std::nullopt;
}

void Body::set_hydration(std::unique_ptr<hydrate::Hydration> hydration) {
    this->hydration = std::move(hydration);
    signal->modified_external();
}

void Body::clear_hydration() {
    hydration->clear();
    signal->modified_external();
}

const std::vector<data::AtomFF>& Body::get_atoms() const {return atoms;}

data::AtomFF& Body::get_atom(unsigned int index) {return atoms[index];}

const data::AtomFF& Body::get_atom(unsigned int index) const {return atoms[index];}

int Body::get_uid() const {return uid;}

std::size_t Body::size_atom() const {return atoms.size();}

std::size_t Body::size_water() const {
    auto h = dynamic_cast<hydrate::ExplicitHydration*>(hydration.get());
    if (h) {return h->waters.size();}
    return 0;
}

std::size_t Body::size_symmetry() const {return symmetries->get().size();}

std::size_t Body::size_symmetry_total() const {
    return std::accumulate(symmetries->get().begin(), symmetries->get().end(), 0, [] (int sum, const symmetry::Symmetry& sym) {return sum + sym.repetitions;});
}