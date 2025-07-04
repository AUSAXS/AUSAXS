// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/generation/JanHydration.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/Molecule.h>
#include <constants/Constants.h>
#include <settings/MoleculeSettings.h>

#include <cassert>

using namespace ausaxs;

hydrate::JanHydration::JanHydration(observer_ptr<data::Molecule> protein) : GridBasedHydration(protein) {
    initialize();
}

hydrate::JanHydration::JanHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy) : GridBasedHydration(protein, std::move(culling_strategy)) {
    initialize();
}

std::span<grid::GridMember<data::Water>> hydrate::JanHydration::generate_explicit_hydration(std::span<grid::GridMember<data::AtomFF>>) {
    assert(protein != nullptr && "JanHydration::generate_explicit_hydration: protein is nullptr.");
    auto grid = protein->get_grid();
    assert(grid != nullptr && "JanHydration::generate_explicit_hydration: grid is nullptr.");

    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    std::size_t water_start = grid->w_members.size();

    auto add_loc = [&] (const Vector3<int>& v) {
        data::Water a(grid->to_xyz(v));
        grid::GridMember<data::Water> gm = grid->add(std::move(a), true);
    };

    // loop over the location of all member atoms
    int r_eff = (grid->get_atomic_radius(form_factor::form_factor_t::C) + grid->get_hydration_radius() + settings::hydrate::shell_correction)/grid->get_width();
    auto[min, max] = grid->bounding_box_index();
    for (int i = min.x(); i < max.x(); i++) {
        int im = std::max(i-r_eff, 0), ip = std::min(i+r_eff, bins.x()-1); // xminus and xplus

        for (int j = min.y(); j < max.y(); j++) {
            int jm = std::max(j-r_eff, 0), jp = std::min(j+r_eff, bins.y()-1); // yminus and yplus

            for (int k = min.z(); k < max.z(); k++) {
                if (gref.is_only_empty_or_volume(i, j, k)) {continue;}
                int km = std::max(k-r_eff, 0), kp = std::min(k+r_eff, bins.z()-1); // zminus and zplus

                // check collisions for x ± r_eff                
                if (gref.is_only_empty_or_volume(im, j, k)) {add_loc(Vector3<int>(im, j, k));}
                if (gref.is_only_empty_or_volume(ip, j, k)) {add_loc(Vector3<int>(ip, j, k));}

                // check collisions for y ± r_eff
                if (gref.is_only_empty_or_volume(i, jp, k)) {add_loc(Vector3<int>(i, jp, k));}
                if (gref.is_only_empty_or_volume(i, jm, k)) {add_loc(Vector3<int>(i, jm, k));}

                // check collisions for z ± r_eff
                if (gref.is_only_empty_or_volume(i, j, km)) {add_loc(Vector3<int>(i, j, km));}
                if (gref.is_only_empty_or_volume(i, j, kp)) {add_loc(Vector3<int>(i, j, kp));}
            }
        }
    }

    return {grid->w_members.begin() + water_start, grid->w_members.end()};
}