#include "PlacementStrategy.h"
#include "Grid.h"

#include "TVector3.h"
#include <iomanip>

/**
 * @brief Description
 */
class RadialPlacement : public PlacementStrategy {
public:
    RadialPlacement(Grid* grid) : PlacementStrategy(grid) {
        prepare_rotations();
    }

    void prepare_rotations(const int divisions = 8) {
        const double w = grid->get_width();
        const int rh = grid->rh, ra = grid->ra;

        vector<vector<int>> bins_1rh;
        vector<vector<int>> bins_2rh;
        vector<vector<int>> bins_rarh;
        double ang = 2*M_PI/divisions;

        for (double theta = 0; theta < 2*M_PI; theta+=ang) {
            for (double phi = 0; phi < 2*M_PI; phi+=ang) {
                double x = cos(phi)*sin(theta);
                double y = sin(phi)*sin(theta);
                double z = cos(theta);
                // cout << (abs(x) < 1e-6 ? 0 : x) << "\t" << (abs(y) < 1e-6 ? 0 : y) << "\t" << (abs(z) < 1e-6 ? 0 : z) << endl;
                // cout << std::setprecision(3) << std::setw(6) << phi/M_PI << ", " << std::setw(6) << theta/M_PI << ": ("
                //      << std::setw(6) << (abs(x) < 1e-6 ? 0 : x) << ", " << std::setw(6) << (abs(y) < 1e-6 ? 0 : y) << ", " << std::setw(6) << (abs(z) < 1e-6 ? 0 : z) << ")" << endl;

                // we use "trunc" for rounding since it is more stable than "round", which often missed some entries due to floating-point errors. 
                bins_1rh.push_back({(int) std::trunc(rh*x/w), (int) std::trunc(rh*y/w), (int) std::trunc(rh*z/w)});
                bins_2rh.push_back({(int) std::trunc(2*(rh-1)*x/w), (int) std::trunc(2*(rh-1)*y/w), (int) std::trunc(2*(rh-1)*z/w)});
                bins_rarh.push_back({(int) std::trunc((ra+rh)*x/w), (int) std::trunc((ra+rh)*y/w), (int) std::trunc((ra+rh)*z/w)});
            }
        }

        // remove duplicate entries
        bins_1rh.erase(std::unique(bins_1rh.begin(), bins_1rh.end()), bins_1rh.end());
        bins_2rh.erase(std::unique(bins_2rh.begin(), bins_2rh.end()), bins_2rh.end());
        bins_rarh.erase(std::unique(bins_rarh.begin(), bins_rarh.end()), bins_rarh.end());
        
        for (const auto& rot : bins_rarh) {
            cout << "" << rot[0] << "\t" << rot[1] << "\t" << rot[2] << "" << endl;
        }

        rot_bins_1rh = bins_1rh;
        rot_bins_2rh = bins_2rh;
        rot_bins_rarh = bins_rarh;
    }

    /**
     * @brief Calculate the bin locations of the rotations.
     * @param ang_divisions the number of angular divisions.
     */
    void prepare_rotations2(const int ang_divisions = 12) {
        const double w = grid->get_width();
        const int rh = grid->rh;

        TVector3 v = TVector3({0, 0, double(rh)});
        vector<vector<int>> bins_1r;
        vector<vector<int>> bins_2r;
        double ang = 2*M_PI/ang_divisions;
        for (int j = 0; j < ang_divisions; j++) {
            if (j == ang_divisions/4) {
                bins_1r.push_back({0, -1*rh, 0});
                bins_2r.push_back({0, -2*rh, 0});
                v.RotateX(ang);
                continue;
            }
            if (j == 3*ang_divisions/4) {
                bins_1r.push_back({0, 1*rh, 0});
                bins_2r.push_back({0, 2*rh, 0});
                v.RotateX(ang);
                continue;
            }
            for (int k = 0; k < ang_divisions; k++) {
                // v.Print();
                bins_1r.push_back({(int) std::lrint(v[0]/w), (int) std::lrint(v[1]/w), (int) std::lrint(v[2]/w)});
                bins_2r.push_back({(int) std::lrint(2*v[0]/w), (int) std::lrint(2*v[1]/w), (int) std::lrint(2*v[2]/w)});
                v.RotateY(ang);
            }
            v.RotateX(ang);
        }

        rot_bins_1rh = bins_1r;
        rot_bins_2rh = bins_2r;

        for (const auto& rot : rot_bins_1rh) {
            cout << "(" << rot[0] << ", " << rot[1] << ", " << rot[2] << ")" << endl;
        }
    }

    vector<vector<int>> rot_bins_1rh; // the bin offsets for rotations of a 1rh length rod
    vector<vector<int>> rot_bins_2rh; // the bin offsets for rotations of a 2rh length rod
    vector<vector<int>> rot_bins_rarh; // the bin offsets for rotations of a 2rh length rod

    vector<shared_ptr<Hetatom>> place() const override {
        // dereference the values we'll need for better performance
        const vector<int> bins = grid->get_bins();
        vector<vector<vector<char>>>& gref = grid->grid;

        // we define a helper lambda
        vector<shared_ptr<Hetatom>> placed_water;
        auto add_loc = [&] (const vector<int> v) {
            shared_ptr<Hetatom> a = Hetatom::create_new_water(grid->to_xyz(v));
            grid->add(a);
            grid->expand_volume(a);
            placed_water.push_back(a);
        };

        for (auto const& pair : grid->members) {
            const vector<int>& loc = pair.second;

            for (auto const& rot : rot_bins_rarh) {
                int xr = loc[0] + rot[0], yr = loc[1] + rot[1], zr = loc[2] + rot[2]; // new coordinates
                
                if (rot == vector({-5, -3, 0})) {
                    cout << "FOUND" << endl;
                }

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

        cout << "TEST:" << endl;
        collision_check({5, 7, 10});
        cout << "GRID: " << gref[5][7][10] << endl;

        return placed_water;
    }

private:
    bool collision_check(const vector<int> loc) const override {
        // dereference the values we'll need for better performance
        vector<vector<vector<char>>>& gref = grid->grid;
        const vector<int> bins = grid->get_bins();
        const int ra = grid->ra, rh = grid->rh;

        // check for collisions at 1rh
        for (auto const& rot : rot_bins_1rh) {
            int xr = loc[0] + rot[0], yr = loc[1] + rot[1], zr = loc[2] + rot[2]; // new coordinates
            if (xr < 0) xr = 0;
            if (xr >= bins[0]) xr = bins[0]-1;
            if (yr < 0) yr = 0;
            if (yr >= bins[1]) yr = bins[1]-1;
            if (zr < 0) zr = 0;
            if (zr >= bins[2]) zr = bins[2]-1;

            if (gref[xr][yr][zr] != 0) {
                // cout << "Collision1 at (" << xr << ", " << yr << ", " << zr << ")" << endl;
                return false;
            };
        }

        // // check for collisions at 2rh
        // for (auto const& rot : rot_bins_2rh) {
        //     int xr = loc[0] + rot[0], yr = loc[1] + rot[1], zr = loc[2] + rot[2]; // new coordinates
        //     if (xr < 0) xr = 0;
        //     if (xr >= bins[0]) xr = bins[0]-1;
        //     if (yr < 0) yr = 0;
        //     if (yr >= bins[1]) yr = bins[1]-1;
        //     if (zr < 0) zr = 0;
        //     if (zr >= bins[2]) zr = bins[2]-1;

        //     if (gref[xr][yr][zr] != 0) {
        //         cout << "Collision2 at (" << xr << ", " << yr << ", " << zr << ")" << endl;
        //         return false;
        //     };
        // }
        // cout << "No collisions!" << endl;
        return true;
    }
};