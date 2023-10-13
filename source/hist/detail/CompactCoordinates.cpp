#include <hist/detail/CompactCoordinates.h>
#include <math/Vector3.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>

using namespace hist::detail;

CompactCoordinates::Data::Data() = default;
CompactCoordinates::Data::Data(const Vector3<double>& v, float w) : x(v.x()), y(v.y()), z(v.z()), w(w) {}

CompactCoordinates::CompactCoordinates(const data::Body& body) : size(body.get_atoms().size()), data(size) {
    for (unsigned int i = 0; i < size; ++i) {
        const auto& a = body.get_atom(i); 
        data[i] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<data::Body>& bodies) {
    size = 0;
    for (const auto& body : bodies) {
        size += body.get_atoms().size();
    }
    data.resize(size);
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i++] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
        }
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<data::record::Water>& atoms) : size(atoms.size()), data(size) {
    for (unsigned int i = 0; i < size; ++i) {
        const auto& a = atoms[i]; 
        data[i] = CompactCoordinates::Data(a.coords, a.effective_charge*a.occupancy);
    }
}

unsigned int CompactCoordinates::get_size() const {return size;}

CompactCoordinates::Data& CompactCoordinates::operator[](unsigned int i) {return data[i];}

const CompactCoordinates::Data& CompactCoordinates::operator[](unsigned int i) const {return data[i];}