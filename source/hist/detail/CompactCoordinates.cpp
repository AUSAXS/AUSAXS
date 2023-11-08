#include <hist/detail/CompactCoordinates.h>
#include <math/Vector3.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>
#include <constants/Constants.h>

using namespace hist::detail;

CompactCoordinates::CompactCoordinates(unsigned int size) : data(size) {}

CompactCoordinates::CompactCoordinates(std::vector<Vector3<double>>&& coordinates, double weight) : data(coordinates.size()) {
    std::transform(coordinates.begin(), coordinates.end(), data.begin(), [weight] (const Vector3<double>& v) {return CompactCoordinatesData(v, weight);});
}

CompactCoordinates::CompactCoordinates(const data::Body& body) : data(body.atom_size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = body.get_atom(i); 
        data[i] = CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<data::Body>& bodies) : data(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.atom_size();})) {
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i++] = CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
        }
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<data::record::Water>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
    }
}

const std::vector<CompactCoordinatesData>& CompactCoordinates::get_data() const {return data;}

std::size_t CompactCoordinates::size() const {return data.size();}

CompactCoordinatesData& CompactCoordinates::operator[](unsigned int i) {return data[i];}

const CompactCoordinatesData& CompactCoordinates::operator[](unsigned int i) const {return data[i];}