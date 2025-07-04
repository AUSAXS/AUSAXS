// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/generation/AxesHydration.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/Molecule.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

#include <cassert>

using namespace ausaxs;

hydrate::AxesHydration::AxesHydration(observer_ptr<data::Molecule> protein) : GridBasedHydration(protein) {
    initialize();
}

hydrate::AxesHydration::AxesHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy) : GridBasedHydration(protein, std::move(culling_strategy)) {
    initialize();
}

hydrate::AxesHydration::~AxesHydration() = default;

void hydrate::AxesHydration::initialize() {
    hydrate::GridBasedHydration::initialize();
}

std::span<grid::GridMember<data::Water>> hydrate::AxesHydration::generate_explicit_hydration(std::span<grid::GridMember<data::AtomFF>> atoms) {
    assert(protein != nullptr && "AxesHydration::generate_explicit_hydration: protein is nullptr.");
    auto grid = protein->get_grid();
    assert(grid != nullptr && "AxesHydration::generate_explicit_hydration: grid is nullptr.");

    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    // short lambda to actually place the generated water molecules
    auto add_loc = [&] (Vector3<double>&& exact_loc) {
        data::Water a(std::move(exact_loc));
        grid::GridMember<data::Water> gm = grid->add(std::move(a), true);
    };

    // loop over the location of all member atoms
    std::size_t water_start = grid->w_members.size();
    double rh = grid->get_hydration_radius() + settings::hydrate::shell_correction;
    for (const auto& atom : atoms) {
        double ra = grid->get_atomic_radius(atom.get_atom_type()); // radius of the atom
        double r_eff_real = ra+rh; // the effective bin radius
        // int r_eff_bin = std::round(r_eff_real)/grid->get_width(); // the effective bin radius in bins

        const auto& coords_abs = atom.get_atom().coordinates();
        auto [x, y, z] = atom.get_bin_loc();

        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
        auto bin_min = grid->to_bins(coords_abs - r_eff_real);
        auto bin_max = grid->to_bins(coords_abs + r_eff_real);
        bin_min.x() = std::max<int>(bin_min.x(), 0); bin_max.x() = std::min<int>(bin_max.x(), bins[0]-1);
        bin_min.y() = std::max<int>(bin_min.y(), 0); bin_max.y() = std::min<int>(bin_max.y(), bins[1]-1);
        bin_min.z() = std::max<int>(bin_min.z(), 0); bin_max.z() = std::min<int>(bin_max.z(), bins[2]-1);

        // check collisions for x ± r_eff
        if ((gref.is_only_empty_or_volume(bin_min.x(), y, z)) && collision_check(Vector3<unsigned int>(bin_min.x(), y, z), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.x() -= r_eff_real;
            add_loc(std::move(exact_loc));
        }
        if ((gref.is_only_empty_or_volume(bin_max.x(), y, z)) && collision_check(Vector3<unsigned int>(bin_max.x(), y, z), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.x() += r_eff_real;
            add_loc(std::move(exact_loc));
        }

        // check collisions for y ± r_eff
        if ((gref.is_only_empty_or_volume(x, bin_min.y(), z)) && collision_check(Vector3<unsigned int>(x, bin_min.y(), z), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.y() -= r_eff_real;
            add_loc(std::move(exact_loc));
        }

        if ((gref.is_only_empty_or_volume(x, bin_max.y(), z)) && collision_check(Vector3<unsigned int>(x, bin_max.y(), z), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.y() += r_eff_real;
            add_loc(std::move(exact_loc));
        }

        // check collisions for z ± r_eff
        if ((gref.is_only_empty_or_volume(x, y, bin_min.z())) && collision_check(Vector3<unsigned int>(x, y, bin_min.z()), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.z() -= r_eff_real;
            add_loc(std::move(exact_loc));
        }

        if ((gref.is_only_empty_or_volume(x, y, bin_max.z())) && collision_check(Vector3<unsigned int>(x, y, bin_max.z()), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.z() += r_eff_real;
            add_loc(std::move(exact_loc));
        }
    }
    return {grid->w_members.begin() + water_start, grid->w_members.end()};
}

bool hydrate::AxesHydration::collision_check(const Vector3<unsigned int>& loc, double ra) const {
    static double rh = constants::radius::get_vdw_radius(constants::atom_t::O); // radius of a water molecule
    auto grid = protein->get_grid();
    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    
    int x = loc.x(), y = loc.y(), z = loc.z();

    // loop over the box [x-r, x+r][y-r, y+r][z-r, z+r]
    int r = gref.is_atom_center(x, y, z)*ra + gref.is_water_center(x, y, z)*rh;

    // we use the range (x-r) to (x+r+1) since the first is inclusive and the second is exclusive. 
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) bins[0])-1; // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) bins[1])-1; // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) bins[2])-1; // zminus and zplus
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                if (!gref.is_only_empty_or_volume(i, j, k) && std::pow(x-i, 2) + std::pow(y-j, 2) + std::pow(z-k, 2) < r*r) {return false;}
            }
        }
    }
    return true;
}