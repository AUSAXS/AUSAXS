#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/detail/CompactCoordinatesData.h>
#include <form_factor/FormFactor.h>
#include <math/Vector3.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>

using namespace hist::detail;

CompactCoordinatesFF::CompactCoordinatesFF(const data::Body& body) : CompactCoordinates(body.atom_size()), ff_types(body.atom_size()) {
    for (unsigned int i = 0; i < size; ++i) {
        const auto& a = body.get_atom(i); 
        data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
        ff_types[i] = static_cast<int>(form_factor::get_type(a.get_element()));
    }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::Body>& bodies) {
    size = std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.atom_size();});
    data.resize(size);
    ff_types.resize(size);
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
            ff_types[i++] = static_cast<int>(form_factor::get_type(a.get_element()));
        }
    }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::record::Water>& atoms) : CompactCoordinates(atoms.size()), ff_types(atoms.size()) {
    for (unsigned int i = 0; i < size; ++i) {
        const auto& a = atoms[i]; 
        data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
        ff_types[i++] = static_cast<int>(form_factor::get_type(a.get_element()));
    }
}

unsigned int CompactCoordinatesFF::get_ff_type(unsigned int i) const {return ff_types[i];}