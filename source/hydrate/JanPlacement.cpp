#include "hydrate/PlacementStrategy.h"
#include "hydrate/Grid.h"

/**
 * @brief JanStrategy
 * 
 * This strategy iterates through all bins, and for every bin which is part of the volume of an atom, it attempts to place a
 * water molecule at x±r, y±r, and z±r. If the location is valid, the molecule will be placed. This will typically generate
 * a lot of molecules, and so a culling method may be useful afterwards. 
 * 
 * It only uses the atomic radius @a ra, and ignores the value set to the hydration atoms @a ra. 
 * 
 */
class JanPlacement : public PlacementStrategy {
public:
    using PlacementStrategy::PlacementStrategy; // inherit constructor
    ~JanPlacement() override {}

    vector<Hetatom> place() const override {
        // dereference the values we'll need for better performance
        vector<vector<vector<char>>>& gref = grid->grid;
        const vector<int> bins = grid->get_bins();

        // place a water molecule (note: not added to the grid before the end of this method)
        vector<Hetatom> placed_water;
        auto add_loc = [&] (const vector<int> v) {
            Hetatom a = Hetatom::create_new_water(grid->to_xyz(v));
            placed_water.push_back(a);
        };

        // loop over the location of all member atoms
        int r_eff = grid->ra;
        vector<vector<int>> bounds = grid->bounding_box();
        for (int i = bounds[0][0]; i < bounds[0][1]; i++) {
            for (int j = bounds[1][0]; j < bounds[1][1]; j++) {
                for (int k = bounds[2][0]; k < bounds[2][1]; k++) {
                    if (gref[i][j][k] == 0) {continue;}

                    // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
                    int im = std::max(i-r_eff, 0), ip = std::min(i+r_eff, bins[0]-1); // xminus and xplus
                    int jm = std::max(j-r_eff, 0), jp = std::min(j+r_eff, bins[1]-1); // yminus and yplus
                    int km = std::max(k-r_eff, 0), kp = std::min(k+r_eff, bins[2]-1); // zminus and zplus

                    // check collisions for x ± r_eff
                    if (gref[im][j][k] == 0) {add_loc({im, j, k});}
                    if (gref[ip][j][k] == 0) {add_loc({ip, j, k});}

                    // check collisions for y ± r_eff
                    if (gref[i][jm][k] == 0) {add_loc({i, jm, k});}
                    if (gref[i][jp][k] == 0) {add_loc({i, jp, k});}

                    // check collisions for z ± r_eff
                    if (gref[i][j][km] == 0) {add_loc({i, j, km});}
                    if (gref[i][j][kp] == 0) {add_loc({i, j, kp});}
                }
            }
        }

        // finally we can add the atoms to the grid
        grid->add(placed_water);
        grid->expand_volume();

        return placed_water;
    }
};