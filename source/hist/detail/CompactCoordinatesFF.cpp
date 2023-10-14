#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactor.h>
#include <math/Vector3.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>

using namespace hist::detail;

CompactCoordinatesFF::Data::Data() = default;
CompactCoordinatesFF::Data::Data(const Vector3<double>& v, float w, form_factor::form_factor_t ff_type) 
    : x(v.x()), y(v.y()), z(v.z()), w(w), ff_type(static_cast<unsigned int>(ff_type)) {}

CompactCoordinatesFF::CompactCoordinatesFF(const data::Body& body) : size(body.get_atoms().size()), data(size) {
    for (unsigned int i = 0; i < size; ++i) {
        const auto& a = body.get_atom(i); 
        data[i] = CompactCoordinatesFF::Data(a.coords, a.effective_charge*a.occupancy, form_factor::get_type(a.get_element()));
    }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::Body>& bodies) {
    size = 0;
    for (const auto& body : bodies) {
        size += body.get_atoms().size();
    }
    data.resize(size);
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i++] = CompactCoordinatesFF::Data(a.coords, a.effective_charge*a.occupancy, form_factor::get_type(a.get_element()));
        }
    }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::record::Water>& atoms) : size(atoms.size()), data(size) {
    for (unsigned int i = 0; i < size; ++i) {
        const auto& a = atoms[i]; 
        data[i] = CompactCoordinatesFF::Data(a.coords, a.effective_charge*a.occupancy, form_factor::get_type(a.get_element()));
    }
}

unsigned int CompactCoordinatesFF::get_size() const {return size;}

CompactCoordinatesFF::Data& CompactCoordinatesFF::operator[](unsigned int i) {return data[i];}

const CompactCoordinatesFF::Data& CompactCoordinatesFF::operator[](unsigned int i) const {return data[i];}