#include "hydrate/PlacementStrategy.h"
#include "hydrate/Grid.h"

/**
 * @brief This strategy iterates through all bins, and for every bin which is part of the volume of an atom, it attempts to place a
 *        water molecule at x±r, y±r, and z±r. If the location is valid, the molecule will be placed. This will typically generate
 *        a lot of molecules, and so a culling method may be useful afterwards. 
 */
class AxesPlacement : public PlacementStrategy {
public:
    using PlacementStrategy::PlacementStrategy; // inherit constructor
    ~AxesPlacement() override {}

    vector<Hetatom> place() const override {
        // dereference the values we'll need for better performance
        vector<vector<vector<char>>>& gref = grid->grid;
        const vector<int> bins = grid->get_bins();
        const int ra = grid->ra; const int rh = grid->rh;

        // we define two helper functions so I can make the checks in the inner loop one-liners
        vector<Hetatom> placed_water;
        auto add_loc = [&] (const vector<int> v) {
            Hetatom a = Hetatom::create_new_water(grid->to_xyz(v));
            grid->add(a);
            grid->expand_volume(a);
            placed_water.push_back(a);
        };

        // loop over the location of all member atoms
        int r_eff = ra+rh;
        vector<Atom> atoms = grid->get_protein_atoms();
        for (auto const& a : atoms) {
            const vector<int>& loc = grid->members.at(a);
            const int x = loc[0], y = loc[1], z = loc[2];

            // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
            int xm = std::max(x-r_eff, 0), xp = std::min(x+r_eff, bins[0]-1); // xminus and xplus
            int ym = std::max(y-r_eff, 0), yp = std::min(y+r_eff, bins[1]-1); // yminus and yplus
            int zm = std::max(z-r_eff, 0), zp = std::min(z+r_eff, bins[2]-1); // zminus and zplus

            // check collisions for x ± r_eff
            if ((gref[xm][y][z] == 0) && collision_check({xm, y, z})) {add_loc({xm, y, z});}
            if ((gref[xp][y][z] == 0) && collision_check({xp, y, z})) {add_loc({xp, y, z});}

            // check collisions for y ± r_eff
            if ((gref[x][ym][z] == 0) && collision_check({x, ym, z})) {add_loc({x, ym, z});}
            if ((gref[x][yp][z] == 0) && collision_check({x, yp, z})) {add_loc({x, yp, z});}

            // check collisions for z ± r_eff
            if ((gref[x][y][zm] == 0) && collision_check({x, y, zm})) {add_loc({x, y, zm});}
            if ((gref[x][y][zp] == 0) && collision_check({x, y, zp})) {add_loc({x, y, zp});}
        }

        return placed_water;
    }

private:
    /**
     * @brief Check if a water molecule can be placed at the given location. 
     * @param loc the location to be checked. 
     * @return True if this is an acceptable location, false otherwise.
     */
    inline bool collision_check(const vector<int> loc) const {
        // dereference the values we'll need for better performance
        vector<vector<vector<char>>>& gref = grid->grid;
        const vector<int> bins = grid->get_bins();
        const int ra = grid->ra, rh = grid->rh;
        const int x = loc[0], y = loc[1], z = loc[2];

        // loop over the box [x-r, x+r][y-r, y+r][z-r, z+r]
        int r = gref[x][y][z] == 'A' ? ra : rh;

        // we use the range (x-r) to (x+r+1) since the first is inclusive and the second is exclusive. 
        int xm = std::max(x-r, 0), xp = std::min(x+r+1, bins[0]-1); // xminus and xplus
        int ym = std::max(y-r, 0), yp = std::min(y+r+1, bins[1]-1); // yminus and yplus
        int zm = std::max(z-r, 0), zp = std::min(z+r+1, bins[2]-1); // zminus and zplus
        for (int i = xm; i < xp; i++) {
            for (int j = ym; j < yp; j++) {
                for (int k = zm; k < zp; k++) {
                    if (gref[i][j][k] != 0 && sqrt(pow(x-i, 2) + pow(y-j, 2) + pow(z-k, 2)) < r) {return false;}
                }
            }
        }
        return true;
    }
};