#include "PlacementStrategy.h"
#include "Grid.h"

/**
 * @brief Description
 */
class RadialPlacement : public PlacementStrategy {
public:
    RadialPlacement(Grid* grid) : PlacementStrategy(grid) {
        prepare_rotations();
    }

    /**
     * @brief Calculate the bin locations of the rotations.
     * @param ang_divisions the number of angular divisions.
     */
    void prepare_rotations(const int ang_divisions = 12) {
        TVector3 v = TVector3({1, 0, 0});
        const double w = grid->get_width();
        const double rh = grid->rh;

        vector<vector<int>> bins_1r;
        vector<vector<int>> bins_2r;
        double ang = 2*M_PI/ang_divisions;
        for (int j = 0; j <= ang_divisions; j++) {
            if (j == ang_divisions/4) {
                bins_1r.push_back({0, 1, 0});
                bins_2r.push_back({0, 2, 0});
                v.RotateY(ang);
                continue;
            }
            if (j == 3*ang_divisions/4) {
                bins_1r.push_back({0, -1, 0});
                bins_2r.push_back({0, -2, 0});
                v.RotateY(ang);
                continue;
            }
            for (int k = 0; k <= ang_divisions; k++) {
                if (k == ang_divisions/4) {
                    bins_1r.push_back({0, 0, 1});
                    bins_2r.push_back({0, 0, 2});
                    v.RotateZ(ang);
                    continue;
                }
                if (k == 3*ang_divisions/4) {
                    bins_1r.push_back({0, 0, -1});
                    bins_2r.push_back({0, 0, -2});
                    v.RotateZ(ang);
                    continue;
                }
                v.Print();
                bins_1r.push_back({(int) std::lrint(rh*v[0]/w), (int) std::lrint(rh*v[1]/w), (int) std::lrint(rh*v[2]/w)});
                bins_2r.push_back({(int) std::lrint(2*rh*v[0]/w), (int) std::lrint(2*rh*v[1]/w), (int) std::lrint(2*rh*v[2]/w)});
                v.RotateZ(ang);
            }
            v.RotateY(ang);
        }

        rot_bins_1r = bins_1r;
        rot_bins_2r = bins_2r;
    }

    vector<vector<int>> rot_bins_1r; // the bin offsets for rotations of a 1rh length rod
    vector<vector<int>> rot_bins_2r; // the bin offsets for rotations of a 2rh length rod

    vector<shared_ptr<Atom>> place(const vector<vector<int>> bounds) override {
        // dereference the values we'll need for better performance
        const vector<int> bins = grid->get_bins();
        vector<vector<vector<char>>>& gref = grid->grid;

        // we define a helper lambda
        vector<shared_ptr<Atom>> placed_water;
        auto add_loc = [&] (const vector<int> v) {
            shared_ptr<Atom> a = Atom::create_new_water(grid->to_xyz(v));
            grid->add(a);
            grid->expand_volume(*a.get());
            placed_water.push_back(a);
        };

        for (auto const& pair : grid->members) {
            const vector<int>& loc = pair.second;

            for (auto const& rot : rot_bins_1r) {
                int xr = loc[0] + rot[0], yr = loc[1] + rot[1], zr = loc[2] + rot[2]; // new coordinates
                
                // check bounds
                if (xr < 0) xr = 0;
                if (xr >= bins[0]) xr = bins[0]-1;
                if (yr < 0) yr = 0;
                if (yr >= bins[1]) yr = bins[1]-1;
                if (zr < 0) zr = 0;
                if (zr >= bins[2]) zr = bins[2]-1;

                if (gref[xr][yr][zr] == 0 && collision_check({xr, yr, zr})) {add_loc({xr, yr, zr});};
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

        // check for collisions at 1rh
        for (auto const& rot : rot_bins_1r) {
            int xr = loc[0] + rot[0], yr = loc[1] + rot[1], zr = loc[2] + rot[2]; // new coordinates
            if (gref[xr][yr][zr]) {return false;};
        }

        // check for collisions at 2rh
        for (auto const& rot : rot_bins_2r) {
            int xr = loc[0] + rot[0], yr = loc[1] + rot[1], zr = loc[2] + rot[2]; // new coordinates
            if (gref[xr][yr][zr]) {return false;};            
        }

        return true;
    }
};