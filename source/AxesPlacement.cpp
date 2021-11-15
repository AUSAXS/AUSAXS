#include "PlacementStrategy.h"
#include "Grid.h"

/**
 * @brief This strategy iterates through all bins, and for every bin which is part of the volume of an atom, it attempts to place a
 *        water molecule at x±r, y±r, and z±r. If the location is valid, the molecule will be placed. This will typically generate
 *        a lot of molecules, and so a culling method may be useful afterwards. 
 */
class AxesPlacement : public PlacementStrategy {
public:
    AxesPlacement(Grid* grid) : PlacementStrategy(grid) {}

    vector<shared_ptr<Atom>> place(const vector<vector<int>> bounds) override {
        // dereference the values we'll need for better performance
        vector<vector<vector<char>>>& gref = grid->grid;
        const vector<int> bins = grid->get_bins();
        const int ra = grid->ra; const int rh = grid->rh;

        // we define two helper functions so I can make the checks in the inner loop one-liners
        vector<shared_ptr<Atom>> placed_water;
        auto add_loc = [&] (const vector<int> v) {
            shared_ptr<Atom> a = Atom::create_new_water(grid->to_xyz(v));
            grid->add(a);
            grid->expand_volume(*a.get());
            placed_water.push_back(a);
        };

        // loop over the minimum bounding box as found above
        for (int i = bounds[0][0]; i < bounds[0][1]; i++) {
            for (int j = bounds[1][0]; j < bounds[1][1]; j++) {
                for (int k = bounds[2][0]; k < bounds[2][1]; k++) {
                    // if this spot is part of an atom
                    if (gref[i][j][k] == 'a') {
                        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
                        int xm = std::max(i-rh, 0), xp = std::min(i+rh, bins[0]); // xminus and xplus
                        int ym = std::max(j-rh, 0), yp = std::min(j+rh, bins[1]); // yminus and yplus
                        int zm = std::max(k-rh, 0), zp = std::min(k+rh, bins[2]); // zminus and zplus

                        // check collisions for x ± r_eff
                        if ((gref[xm][j][k] == 0) && collision_check({xm, j, k})) add_loc({xm, j, k});
                        if ((gref[xp][j][k] == 0) && collision_check({xp, j, k})) add_loc({xp, j, k});

                        // check collisions for y ± r_eff
                        if ((gref[i][ym][k] == 0) && collision_check({i, ym, k})) add_loc({i, ym, k});
                        if ((gref[i][yp][k] == 0) && collision_check({i, yp, k})) add_loc({i, yp, k});

                        // check collisions for z ± r_eff
                        if ((gref[i][j][zm] == 0) && collision_check({i, j, zm})) add_loc({i, j, zm});
                        if ((gref[i][j][zp] == 0) && collision_check({i, j, zp})) add_loc({i, j, zp});
                    } 
                }
            }
        }
        return placed_water;
    }

private:
    bool collision_check(const vector<int> loc) const override {
        // dereference the values we'll need for better performance
        vector<vector<vector<char>>>& gref = grid->grid;
        const vector<int> bins = grid->get_bins();
        const int ra = grid->ra, rh = grid->rh;
        const int x = loc[0], y = loc[1], z = loc[2];

        // loop over the box [x-r, x+r][y-r, y+r][z-r, z+r]
        int r = gref[x][y][z] == 'A' ? ra : rh;

        // we use the range (x-r+1) to (x+r) since the first is inclusive and the second is exclusive. 
        // thus the range that is being checked is actually [x-r+1, x+r-1] because we allow the surfaces to overlap
        int xm = std::max(x-r+1, 0), xp = std::min(x+r, bins[0]); // xminus and xplus
        int ym = std::max(y-r+1, 0), yp = std::min(y+r, bins[1]); // yminus and yplus
        int zm = std::max(z-r+1, 0), zp = std::min(z+r, bins[2]); // zminus and zplus
        for (int i = xm; i < xp; i++) {
            for (int j = ym; j < yp; j++) {
                for (int k = zm; k < zp; k++) {
                    if (gref[i][j][k] != 0) {return false;}
                }
            }
        }
        return true;
    }
};