// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/detail/SimpleExvHelpers.h>

#include <data/Molecule.h>
#include <data/symmetry/MoleculeSymmetryFacade.h>
#include <form_factor/FormFactorInfo.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/detail/CompactCoordinates.h>
#include <settings/InternalState.h>

using namespace ausaxs;
using namespace ausaxs::hist::detail;

double simple_exv::get_volume_per_atom(observer_ptr<const data::Molecule> molecule) {
    assert(molecule != nullptr && "simple_exv::get_volume_per_atom: molecule is nullptr.");
    int n_atoms = molecule->symmetry().size_atom_total();
    assert(0 < n_atoms && "simple_exv::get_volume_per_atom: Division by zero. The molecule has no atoms.");
    return molecule->get_volume_grid()/static_cast<double>(n_atoms);
}

void simple_exv::apply(hist::detail::CompactCoordinates& data_a, double volume_per_atom) {
    data_a.implicit_excluded_volume(volume_per_atom);
}

std::vector<double> simple_exv::get_type_weights(observer_ptr<const data::Molecule> molecule) {
    double displaced_charge = settings::internal_state::solvent_density()*get_volume_per_atom(molecule);
    const auto& ff_indices = form_factor::manager::get_active_indices();
    int n_ff = form_factor::get_active_count();

    std::vector<double> w(n_ff, 0); // no scatterer has the excluded volume type
    for (int ff = form_factor::start_index_for_explicit_exv(); ff < n_ff; ++ff) {
        w[ff] = constants::charge::get_ff_charge(static_cast<form_factor::form_factor_t>(ff_indices[ff]));
        if (ff != form_factor::water_bin) {w[ff] -= displaced_charge;}
    }
    return w;
}
