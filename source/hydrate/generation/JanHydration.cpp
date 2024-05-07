/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/generation/JanHydration.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/record/Water.h>
#include <constants/Constants.h>

using namespace data::record;

std::vector<data::record::Water> hydrate::JanHydration::generate_explicit_hydration() {
    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    std::vector<Water> placed_water;
    placed_water.reserve(grid->a_members.size());
    auto add_loc = [&] (const Vector3<int>& v) {
        placed_water.emplace_back(Water::create_new_water(grid->to_xyz(v)));
    };

    // loop over the location of all member atoms
    int r_eff = (grid->get_atomic_radius(constants::atom_t::C) + grid->get_hydration_radius())/grid->get_width();
    auto[min, max] = grid->bounding_box_index();
    for (int i = min.x(); i < max.x(); i++) {
        for (int j = min.y(); j < max.y(); j++) {
            for (int k = min.z(); k < max.z(); k++) {
                if (gref.index(i, j, k) == grid::detail::EMPTY) {continue;}

                // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
                int im = std::max(i-r_eff, 0), ip = std::min(i+r_eff, (int) bins.x()-1); // xminus and xplus
                int jm = std::max(j-r_eff, 0), jp = std::min(j+r_eff, (int) bins.y()-1); // yminus and yplus
                int km = std::max(k-r_eff, 0), kp = std::min(k+r_eff, (int) bins.z()-1); // zminus and zplus

                // check collisions for x ± r_eff                
                if (gref.index(im, j, k) == grid::detail::EMPTY) {add_loc(Vector3<int>(im, j, k));}
                if (gref.index(ip, j, k) == grid::detail::EMPTY) {add_loc(Vector3<int>(ip, j, k));}

                // check collisions for y ± r_eff
                if (gref.index(i, jp, k) == grid::detail::EMPTY) {add_loc(Vector3<int>(i, jp, k));}
                if (gref.index(i, jm, k) == grid::detail::EMPTY) {add_loc(Vector3<int>(i, jm, k));}

                // check collisions for z ± r_eff
                if (gref.index(i, j, km) == grid::detail::EMPTY) {add_loc(Vector3<int>(i, j, km));}
                if (gref.index(i, j, kp) == grid::detail::EMPTY) {add_loc(Vector3<int>(i, j, kp));}
            }
        }
    }
    return placed_water;
}