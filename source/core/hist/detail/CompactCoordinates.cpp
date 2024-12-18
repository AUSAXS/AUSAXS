/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/detail/CompactCoordinates.h>
#include <data/Body.h>
#include <constants/Constants.h>

using namespace ausaxs::hist::detail;

CompactCoordinates::CompactCoordinates(unsigned int size) : data(size) {}

CompactCoordinates::CompactCoordinates(std::vector<Vector3<double>>&& coordinates, double weight) : data(coordinates.size()) {
    std::transform(coordinates.begin(), coordinates.end(), data.begin(), [weight] (const Vector3<double>& v) {return CompactCoordinatesData(v, weight);});
}

CompactCoordinates::CompactCoordinates(const std::vector<data::AtomFF>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = CompactCoordinatesData(a.coordinates(), a.weight());
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<data::Body>& bodies) : data(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.size_atom();})) {
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i++] = CompactCoordinatesData(a.coordinates(), a.weight());
        }
    }
}

CompactCoordinates::CompactCoordinates(const std::vector<data::Water>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = CompactCoordinatesData(a.coordinates(), a.weight());
    }
}

void CompactCoordinates::implicit_excluded_volume(double volume_per_atom) {
    double displaced_charge = constants::charge::density::water*volume_per_atom;
    double charge_per_atom = -displaced_charge;
    std::for_each(data.begin(), data.end(), [charge_per_atom] (CompactCoordinatesData& d) {d.value.w += charge_per_atom;});
}

std::vector<CompactCoordinatesData>& CompactCoordinates::get_data() {return data;}

const std::vector<CompactCoordinatesData>& CompactCoordinates::get_data() const {return data;}

std::size_t CompactCoordinates::size() const {return data.size();}

CompactCoordinatesData& CompactCoordinates::operator[](unsigned int i) {return data[i];}

const CompactCoordinatesData& CompactCoordinates::operator[](unsigned int i) const {return data[i];}