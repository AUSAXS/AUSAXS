#include <hydrate/OutlierCulling.h>
#include <hydrate/Grid.h>

#include <utility>

vector<Hetatom> grid::OutlierCulling::cull(vector<grid::GridMember<Hetatom>>& placed_water) const {
    if (target_count == 0) {
        vector<Hetatom> final_water(placed_water.size());
        std::transform(placed_water.begin(), placed_water.end(), final_water.begin(), [] (GridMember<Hetatom>& gm) {return gm.atom;});
        return final_water;
    }

    vector<std::pair<GridMember<Hetatom>, int>> v(placed_water.size());
    const int r = 3*grid->ra; // use 2*atomic_radius as the boundary
    auto bins = grid->get_bins();
    const vector<vector<vector<char>>>& gref = grid->grid;
    size_t index = 0;
    for (const auto& water : placed_water) {
        const int x = water.loc[0], y = water.loc[1], z = water.loc[2];
        int score = 0;

        // create a box of size [x-2r, x+2r][y-2r, y+2r][z-2r, z+2r] within the bounds
        int xm = std::max(x-r, 0), xp = std::min(x+r+1, int(bins[0])-1); // xminus and xplus
        int ym = std::max(y-r, 0), yp = std::min(y+r+1, int(bins[1])-1); // yminus and yplus
        int zm = std::max(z-r, 0), zp = std::min(z+r+1, int(bins[2])-1); // zminus and zplus

        for (int i = xm; i < xp; i++) {
            for (int j = ym; j < yp; j++) {
                for (int k = zm; k < zp; k++) {
                    if (gref[i][j][k] == 'A') {score++;}
                    else if (gref[i][j][k] == 'H') {score-=5;}
                }
            }
        }
        v[index++] = std::make_pair(water, score);
    }

    // sort the scores
    std::sort(v.begin(), v.end(), [](auto &left, auto &right) {return left.second < right.second;});

    // copy the first target_count entries in the sorted vector
    vector<Hetatom> final_water(target_count); // the final water molecules that will be used
    vector<Hetatom> removed_water(placed_water.size() - target_count); // the water molecules which will be removed
    size_t n = 0;
    for (n = 0; n < target_count; n++) {
        final_water[n] = v[n].first.atom;
    }
    for (; n < placed_water.size(); n++) {
        removed_water[n] = v[n].first.atom;
    }

    grid->remove(removed_water);
    return final_water;
}