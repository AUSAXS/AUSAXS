#include <hist/detail/CompactCoordinates.h>
#include <math/Vector3.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <data/Body.h>

using namespace hist::detail;

CompactCoordinates::Data::Data() = default;
CompactCoordinates::Data::Data(const Vector3<double>& v, float w) : x(v.x()), y(v.y()), z(v.z()), w(w) {}

CompactCoordinates::CompactCoordinates(const Body& body) : size(body.get_atoms().size()), data(size) {
    for (unsigned int i = 0; i < size; ++i) {
        const Atom& a = body.get_atom(i); 
        data[i] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<Body>& bodies) {
    size = 0;
    for (const Body& body : bodies) {
        size += body.get_atoms().size();
    }
    data.resize(size);
    unsigned int i = 0;
    for (const Body& body : bodies) {
        for (const Atom& a : body.get_atoms()) {
            data[i++] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
        }
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<Water>& atoms) : size(atoms.size()), data(size) {
    for (unsigned int i = 0; i < size; ++i) {
        const Water& a = atoms[i]; 
        data[i] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
    }
}

unsigned int CompactCoordinates::get_size() const {return size;}

CompactCoordinates::Data& CompactCoordinates::operator[](unsigned int i) {return data[i];}

const CompactCoordinates::Data& CompactCoordinates::operator[](unsigned int i) const {return data[i];}