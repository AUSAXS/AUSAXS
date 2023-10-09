#include <hydrate/culling/OutlierCulling.h>
#include <hydrate/GridMember.h>
#include <hydrate/Grid.h>
#include <math/Vector3.h>
#include <data/Water.h>
#include <utility/Constants.h>

#include <utility>

std::vector<Water> grid::OutlierCulling::cull(std::vector<grid::GridMember<Water>>& placed_water) const {
    if (target_count == 0) {
        std::vector<Water> final_water(placed_water.size());
        std::transform(placed_water.begin(), placed_water.end(), final_water.begin(), [] (GridMember<Water>& gm) {return gm.get_atom();});
        return final_water;
    }

    std::vector<std::pair<GridMember<Water>, int>> v(placed_water.size());
    int r = 3*grid->get_atomic_radius(constants::atom_t::C)/grid->get_width(); // use 2*atomic_radius as the boundary
    auto bins = grid->get_bins();
    const GridObj& gref = grid->grid;
    unsigned int index = 0;
    for (const auto& water : placed_water) {
        const int x = water.get_loc().x(), y = water.get_loc().y(), z = water.get_loc().z();
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
    std::vector<Water> final_water(target_count); // the final water molecules that will be used
    std::vector<Water> removed_water(placed_water.size() - target_count); // the water molecules which will be removed
    unsigned int n = 0;
    for (n = 0; n < target_count; n++) {
        final_water[n] = v[n].first.get_atom();
    }
    for (; n < placed_water.size(); n++) {
        removed_water[n] = v[n].first.get_atom();
    }

    grid->remove(removed_water);
    return final_water;
}