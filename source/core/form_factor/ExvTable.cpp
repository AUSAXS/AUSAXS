// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/ExvTable.h>
#include <settings/MoleculeSettings.h>

#include <stdexcept>
#include <string>

using namespace ausaxs;

constants::exv::detail::ExvSet constants::exv::get_exv_set() {
    switch (settings::molecule::exv_set) {
        case settings::molecule::ExvSet::Traube: return Traube;
        case settings::molecule::ExvSet::Voronoi_explicit_H: return Voronoi_explicit_H;
        case settings::molecule::ExvSet::Voronoi_implicit_H: return Voronoi_implicit_H;
        case settings::molecule::ExvSet::MinimumFluctutation_explicit_H: return MinimumFluctuation_explicit_H;
        case settings::molecule::ExvSet::MinimumFluctutation_implicit_H: return MinimumFluctuation_implicit_H;
        case settings::molecule::ExvSet::vdw: return vdw;
        default: 
            throw std::runtime_error(
                "constants::displaced_volume::get_exv_set: Invalid displaced volume set" 
                "(enum " + std::to_string(static_cast<int>(settings::molecule::exv_set)) + ")");
    }
}