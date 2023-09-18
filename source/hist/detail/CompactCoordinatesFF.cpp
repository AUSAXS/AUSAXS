#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/detail/FormFactor.h>
#include <math/Vector3.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <data/Body.h>

using namespace hist::detail;

CompactCoordinatesFF::Data::Data() = default;
CompactCoordinatesFF::Data::Data(const Vector3<double>& v, float w, form_factor_t ff_type) 
    : x(v.x()), y(v.y()), z(v.z()), w(w), ff_type(static_cast<unsigned int>(ff_type)) {}

CompactCoordinatesFF::CompactCoordinatesFF(const Body& body) : size(body.get_atoms().size()), data(size) {
    for (unsigned int i = 0; i < size; ++i) {
        const Atom& a = body.get_atom(i); 
        data[i] = CompactCoordinatesFF::Data(a.coords, a.effective_charge*a.occupancy, FormFactor::get_type(a));
    }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<Body>& bodies) {
    size = 0;
    for (const Body& body : bodies) {
        size += body.get_atoms().size();
    }
    data.resize(size);
    unsigned int i = 0;
    for (const Body& body : bodies) {
        for (const Atom& a : body.get_atoms()) {
            data[i++] = CompactCoordinatesFF::Data(a.coords, a.effective_charge*a.occupancy, FormFactor::get_type(a));
        }
    }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<Water>& atoms) : size(atoms.size()), data(size) {
    for (unsigned int i = 0; i < size; ++i) {
        const Water& a = atoms[i]; 
        data[i] = CompactCoordinatesFF::Data(a.coords, a.effective_charge*a.occupancy, FormFactor::get_type(a));
    }
}