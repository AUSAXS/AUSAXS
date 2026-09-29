// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/CompactCoordinatesFactory.h>

#include <data/Molecule.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::hist::detail;

namespace {
    template<typename Atoms>
    std::vector<CompactCoordinates> construct_by_form_factor(const Atoms& atoms) {
        auto map = form_factor::manager::get_active_mapping();
        auto active_index = [&map] (const data::AtomFF& atom) {
            if (atom.form_factor_type() == form_factor::form_factor_t::UNKNOWN) {
                throw except::runtime_error(
                    "factory::construct<true>: Attempted to use an atom with UNKNOWN form factor type.\n"
                    "Form factor information is required for the selected excluded volume model."
                );
            }
            return map[static_cast<int>(atom.form_factor_type())];
        };

        int form_factor_count = form_factor::get_active_count();
        std::vector<int> counts(form_factor_count, 0);
        for (const auto& atom : atoms) {++counts[active_index(atom)];}

        std::vector<CompactCoordinates> coordinates(form_factor_count);
        for (int form_factor = 0; form_factor < form_factor_count; ++form_factor) {
            coordinates[form_factor].resize(counts[form_factor]);
        }

        std::vector<int> filled(form_factor_count, 0);
        for (const auto& atom : atoms) {
            int form_factor = active_index(atom);
            int index = filled[form_factor]++;
            coordinates[form_factor].set_position(index, atom.coordinates());
            coordinates[form_factor].get_weight(index) = static_cast<float>(atom.weight());
        }
        return coordinates;
    }
}

template<bool form_factors>
factory::AtomicCoordinates<form_factors> factory::construct(const std::vector<data::AtomFF>& atoms) {
    if constexpr (form_factors) {
        return construct_by_form_factor(atoms);
    } else {
        CompactCoordinates coordinates;
        coordinates.fill(atoms);
        return coordinates;
    }
}

template<bool form_factors>
factory::AtomicCoordinates<form_factors> factory::construct_from_atoms(observer_ptr<const data::Molecule> molecule) {
    if constexpr (form_factors) {
        return construct_by_form_factor(molecule->iterate_atoms());
    } else {
        CompactCoordinates coordinates;
        coordinates.fill_from_atoms(molecule);
        return coordinates;
    }
}

template CompactCoordinates factory::construct<false>(const std::vector<data::AtomFF>&);
template std::vector<CompactCoordinates> factory::construct<true>(const std::vector<data::AtomFF>&);
template CompactCoordinates factory::construct_from_atoms<false>(observer_ptr<const data::Molecule>);
template std::vector<CompactCoordinates> factory::construct_from_atoms<true>(observer_ptr<const data::Molecule>);