#include "hydrate/CullingStrategy.h"
#include "hydrate/Grid.h"

#include <utility>

/**
 * @brief Iterate through all water molecules, and count how many other molecules are nearby. Atoms counts as +1, while other water molecules counts as -2. 
 *        Then start removing the most negative water molecules until the desired count is reached. 
 */
class OutlierCulling : public CullingStrategy {
public:
    using CullingStrategy::CullingStrategy;
    ~OutlierCulling() override {}

    // runs in O(n ln n) where n is the number of water molecules
    vector<Hetatom> cull(vector<Hetatom>& placed_water) const override {
        if (target_count == 0) {
            return placed_water;
        }

        vector<std::pair<Hetatom, int>> v(placed_water.size());
        const int r = 3*grid->ra; // use 2*atomic_radius as the boundary
        const vector<int> bins = grid->get_bins();
        const vector<vector<vector<char>>>& gref = grid->grid;
        for (size_t n = 0; n < placed_water.size(); n++) {
            const vector<int> loc = grid->w_members.at(placed_water[n]);
            const int x = loc[0], y = loc[1], z = loc[2];
            int score = 0;

            // create a box of size [x-2r, x+2r][y-2r, y+2r][z-2r, z+2r] within the bounds
            int xm = std::max(x-r, 0), xp = std::min(x+r+1, bins[0]-1); // xminus and xplus
            int ym = std::max(y-r, 0), yp = std::min(y+r+1, bins[1]-1); // yminus and yplus
            int zm = std::max(z-r, 0), zp = std::min(z+r+1, bins[2]-1); // zminus and zplus

            for (int i = xm; i < xp; i++) {
                for (int j = ym; j < yp; j++) {
                    for (int k = zm; k < zp; k++) {
                        if (gref[i][j][k] == 'A') {score++;}
                        else if (gref[i][j][k] == 'H') {score-=5;}
                    }
                }
            }
            v[n] = std::make_pair(placed_water[n], score);
        }

        // sort the scores
        std::sort(v.begin(), v.end(), [](auto &left, auto &right) {return left.second < right.second;});

        // copy the first target_count entries in the sorted vector
        vector<Hetatom> final_water(target_count); // the final water molecules that will be used
        for (int n = 0; n < target_count; n++) {
            final_water[n] = v[n].first;
        }

        // remove the remaining entries
        for (size_t n = target_count; n < v.size(); n++) {
            grid->remove(v[n].first);
        }
        return final_water;
    }
};