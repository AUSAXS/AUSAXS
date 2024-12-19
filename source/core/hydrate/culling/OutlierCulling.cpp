/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "form_factor/FormFactorType.h"
#include <hydrate/culling/OutlierCulling.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/Molecule.h>
#include <constants/Constants.h>

#include <utility>

using namespace ausaxs;
using namespace ausaxs::hydrate;

void OutlierCulling::cull(std::vector<grid::GridMember<data::Water>>& placed_water) const {
    auto grid = molecule->get_grid();

    if (target_count == 0) {return;}

    std::vector<std::pair<grid::GridMember<data::Water>, int>> v(placed_water.size());
    int r = 3*grid->get_atomic_radius(form_factor::form_factor_t::C)/grid->get_width(); // use 2*atomic_radius as the boundary
    auto bins = grid->get_bins();
    const grid::detail::GridObj& gref = grid->grid;
    unsigned int index = 0;
    for (const auto& water : placed_water) {
        const int x = water.get_bin_loc().x(), y = water.get_bin_loc().y(), z = water.get_bin_loc().z();
        int score = 0;

        // create a box of size [x-2r, x+2r][y-2r, y+2r][z-2r, z+2r] within the bounds
        int xm = std::max(x-r, 0), xp = std::min(x+r+1, int(bins[0])-1); // xminus and xplus
        int ym = std::max(y-r, 0), yp = std::min(y+r+1, int(bins[1])-1); // yminus and yplus
        int zm = std::max(z-r, 0), zp = std::min(z+r+1, int(bins[2])-1); // zminus and zplus

        for (int i = xm; i < xp; i++) {
            for (int j = ym; j < yp; j++) {
                for (int k = zm; k < zp; k++) {
                    if (gref.is_atom_center(i, j, k)) {score++;}
                    else if (gref.is_water_center(i, j, k)) {score-=5;}
                }
            }
        }
        v[index++] = std::make_pair(water, score);
    }

    // sort the scores
    std::sort(v.begin(), v.end(), [](auto &left, auto &right) {return left.second < right.second;});

    // copy the first target_count entries in the sorted vector
    std::vector<bool> to_remove(placed_water.size() - target_count, false);
    unsigned int n = 0;
    while (n < target_count) {n++;}
    for (; n < placed_water.size(); n++) {
        to_remove[n] = true;
    }
    grid->remove_waters(to_remove);
}