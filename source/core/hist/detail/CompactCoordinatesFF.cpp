/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/detail/CompactCoordinatesData.h>
#include <form_factor/FormFactor.h>
#include <math/Vector3.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>

using namespace hist::detail;

CompactCoordinatesFF::CompactCoordinatesFF(const data::Body& body) : CompactCoordinates(body.size_atom()), ff_types(body.size_atom()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = body.get_atom(i); 
        data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
        ff_types[i] = static_cast<int>(form_factor::get_type(a.get_element(), a.get_atomic_group()));
    }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::Body>& bodies) : CompactCoordinates(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.size_atom();})) {
    ff_types.resize(size());
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
            ff_types[i++] = static_cast<int>(form_factor::get_type(a.get_element(), a.get_atomic_group()));
        }
    }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::record::Water>& atoms) : CompactCoordinates(atoms.size()), ff_types(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
        ff_types[i] = static_cast<int>(form_factor::get_type(a.get_element(), a.get_atomic_group()));
    }
}

unsigned int CompactCoordinatesFF::get_ff_type(unsigned int i) const {return ff_types[i];}