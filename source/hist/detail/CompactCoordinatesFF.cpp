#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/detail/CompactCoordinatesData.h>
#include <form_factor/FormFactor.h>
#include <math/Vector3.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>

using namespace hist::detail;

CompactCoordinatesFF::CompactCoordinatesFF(const data::Body& body) : CompactCoordinates(body.atom_size()), ff_types(body.atom_size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = body.get_atom(i); 
        data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
        ff_types[i] = static_cast<int>(form_factor::get_type(a.get_element(), a.get_atomic_group()));
    }
}

#include <form_factor/ExvFormFactor.h>
#include <hist/foxs/FormFactorFoXS.h>
#include <fstream>
CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::Body>& bodies) : CompactCoordinates(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.atom_size();})) {
    ff_types.resize(size());
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
            ff_types[i++] = static_cast<int>(form_factor::get_type(a.get_element(), a.get_atomic_group()));
            // ff_types[i++] = static_cast<int>(form_factor::get_type(a.get_element()));

            // auto name = a.get_group_name(); auto residue = a.get_residue_name();
            // if (name == "N") {
            //     if (residue == "PRO") {ff_types[i-1] = static_cast<int>(form_factor::form_factor_t::N);}
            //     else {ff_types[i-1] = static_cast<int>(form_factor::form_factor_t::NH);}
            // }
            // else if (name == "NE2" && residue != "GLN") {ff_types[i-1] = static_cast<int>(form_factor::form_factor_t::N);}
            // else if (name == "OD2") {ff_types[i-1] = static_cast<int>(form_factor::form_factor_t::O);}
            // else if (name == "OE2") {ff_types[i-1] = static_cast<int>(form_factor::form_factor_t::O);}
            // else if (name == "CH2") {ff_types[i-1] = static_cast<int>(form_factor::form_factor_t::CH2);}
        }
    }

    // std::ofstream out("ff.txt");
    // for (unsigned int i = 0; i < 2367; ++i) {
    //     auto fft = static_cast<form_factor::form_factor_t>(ff_types[i]);
    //     out << "i = " << i+1 << ": " << bodies[0].get_atom(i).get_group_name() << " " << form_factor::to_string(fft) << " " << form_factor::foxs::storage::atomic::get_form_factor(fft).evaluate(0) << " " << form_factor::foxs::storage::exv::get_form_factor(fft).evaluate(0) << std::endl;
    // }
}

CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::record::Water>& atoms) : CompactCoordinates(atoms.size()), ff_types(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = hist::detail::CompactCoordinatesData(a.coords, a.effective_charge*a.occupancy);
        ff_types[i] = static_cast<int>(form_factor::get_type(a.get_element(), a.get_atomic_group()));
    }
}

unsigned int CompactCoordinatesFF::get_ff_type(unsigned int i) const {return ff_types[i];}