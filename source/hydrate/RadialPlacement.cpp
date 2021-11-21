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

    // void prepare_rotations(const int divisions = 8) {
    //     const double w = grid->get_width();
    //     const int rh = grid->rh, ra = grid->ra;

    //     // vector<vector<int>> bins;
    //     vector<vector<int>> bins_1rh;
    //     vector<vector<int>> bins_3rh;
    //     vector<vector<int>> bins_5rh;
    //     vector<vector<int>> bins_7rh;
    //     vector<vector<int>> bins_rarh;
    //     double ang = 2*M_PI/divisions;

    //     for (double theta = 0; theta < 2*M_PI; theta+=ang) {
    //         for (double phi = 0; phi < M_PI; phi+=ang) {
    //             double x = cos(phi)*sin(theta);
    //             double y = sin(phi)*sin(theta);
    //             double z = cos(theta);
    //             // cout << (abs(x) < 1e-6 ? 0 : x) << "\t" << (abs(y) < 1e-6 ? 0 : y) << "\t" << (abs(z) < 1e-6 ? 0 : z) << endl;
    //             // cout << std::setprecision(3) << std::setw(6) << phi/M_PI << ", " << std::setw(6) << theta/M_PI << ": ("
    //             //      << std::setw(6) << (abs(x) < 1e-6 ? 0 : x) << ", " << std::setw(6) << (abs(y) < 1e-6 ? 0 : y) << ", " << std::setw(6) << (abs(z) < 1e-6 ? 0 : z) << ")" << endl;

    //             // we use "trunc" for rounding since it is more stable than "round", which often missed some entries due to floating-point errors. 
    //             // bins.push_back({(int) std::trunc(x/w), (int) std::trunc(y/w), (int) std::trunc(z/w)});
    //             bins_1rh.push_back({(int) std::trunc(rh*x/w), (int) std::trunc(rh*y/w), (int) std::trunc(rh*z/w)});
    //             bins_3rh.push_back({(int) std::trunc(3*rh*x/w), (int) std::trunc(3*rh*y/w), (int) std::trunc(3*rh*z/w)});
    //             bins_5rh.push_back({(int) std::trunc(5*rh*x/w), (int) std::trunc(5*rh*y/w), (int) std::trunc(5*rh*z/w)});
    //             bins_7rh.push_back({(int) std::trunc(7*rh*x/w), (int) std::trunc(7*rh*y/w), (int) std::trunc(7*rh*z/w)});
    //             bins_rarh.push_back({(int) std::trunc((ra+rh)*x/w), (int) std::trunc((ra+rh)*y/w), (int) std::trunc((ra+rh)*z/w)});
    //         }
    //     }

    //     // remove duplicate entries
    //     // bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
    //     bins_1rh.erase(std::unique(bins_1rh.begin(), bins_1rh.end()), bins_1rh.end());
    //     bins_3rh.erase(std::unique(bins_3rh.begin(), bins_3rh.end()), bins_3rh.end());
    //     bins_5rh.erase(std::unique(bins_5rh.begin(), bins_5rh.end()), bins_5rh.end());
    //     bins_7rh.erase(std::unique(bins_7rh.begin(), bins_7rh.end()), bins_7rh.end());
    //     bins_rarh.erase(std::unique(bins_rarh.begin(), bins_rarh.end()), bins_rarh.end());
        
    //     for (const auto& rot : bins_rarh) {
    //         cout << "" << rot[0] << "\t" << rot[1] << "\t" << rot[2] << "" << endl;
    //     }

    //     // rot_bins = bins;
    //     rot_bins_1rh = bins_1rh;
    //     rot_bins_3rh = bins_3rh;
    //     rot_bins_5rh = bins_5rh;
    //     rot_bins_7rh = bins_7rh;
    //     rot_bins_rarh = bins_rarh;

    //     // sanity check
    //     if (bins_1rh.size() != bins_3rh.size() || bins_1rh.size() != bins_5rh.size() 
    //         || bins_1rh.size() != bins_7rh.size() || bins_1rh.size() != bins_rarh.size()) {                
    //         print_err("Error in RadialPlacement::prepare_rotations: Sizes of rotation vectors are not equal!");
    //         exit(1);
    //     }
    // }

    void prepare_rotations(const int divisions = 8) {
        const double w = grid->get_width();
        const int rh = grid->rh, ra = grid->ra;

        // vector<vector<int>> bins;
        vector<vector<int>> bins_1rh;
        vector<vector<int>> bins_3rh;
        vector<vector<int>> bins_5rh;
        vector<vector<int>> bins_7rh;
        vector<vector<int>> bins_rarh;
        double ang = 2*M_PI/divisions;

        // we generate one octant of a sphere, and then reflect it to generate the rest
        // we do this to ensure the sphere is symmetric. If we simply generate it all at once, floating-point errors moves some of the bins around
        vector<vector<double>> sphere;
        for (double theta = 0; theta <= M_PI*0.5; theta+=ang) {
            for (double phi = 0; phi <= M_PI*0.5; phi+=ang) {
                double x = cos(phi)*sin(theta);
                double y = sin(phi)*sin(theta);
                double z = cos(theta);
                sphere.push_back({x, y, z});
                sphere.push_back({-x, y, z});
                sphere.push_back({x, -y, z});
                sphere.push_back({-x, -y, z});
                sphere.push_back({x, y, -z});
                sphere.push_back({-x, y, -z});
                sphere.push_back({x, -y, -z});
                sphere.push_back({-x, -y, -z});
            }
        }

        // remove duplicates
        vector<vector<double>> rots;
        for (auto& p : sphere) {
            bool present = false;
            for (int i = 0; i < 3; i++) { // fix the easy floating point errors
                if (abs(p[i]) < 1e-5) {p[i] = 0;}
            }
            for (const auto& r : rots) { // go through all rotations and try to find a duplicate entry
                double sum = 0;
                for (int i = 0; i < 3; i++) { // sum the difference for the vector
                    sum += abs(p[i]-r[i]);
                }
                if (sum < 1e-5) { // if the total difference is small
                    present = true; // the element is already present
                }
            }
            if (!present) { // if the element was not already present
                rots.push_back(p); // add it
            }
        }

        for (const auto& rot : rots) {
            const double xr = rot[0], yr = rot[1], zr = rot[2];
            bins_1rh.push_back({(int) std::trunc(rh*xr), (int) std::trunc(rh*yr), (int) std::trunc(rh*zr)});
            bins_3rh.push_back({(int) std::trunc(3*rh*xr), (int) std::trunc(3*rh*yr), (int) std::trunc(3*rh*zr)});
            bins_5rh.push_back({(int) std::trunc(5*rh*xr), (int) std::trunc(5*rh*yr), (int) std::trunc(5*rh*zr)});
            bins_7rh.push_back({(int) std::trunc(7*rh*xr), (int) std::trunc(7*rh*yr), (int) std::trunc(7*rh*zr)});
            bins_rarh.push_back({(int) std::trunc((ra+rh)*xr), (int) std::trunc((ra+rh)*yr), (int) std::trunc((ra+rh)*zr)});
        }

        // rot_bins = bins;
        rot_bins_1rh = bins_1rh;
        rot_bins_3rh = bins_3rh;
        rot_bins_5rh = bins_5rh;
        rot_bins_7rh = bins_7rh;
        rot_bins_rarh = bins_rarh;
    }

    // vector<vector<int>> rot_bins; // the bin offsets representing rotations
    vector<vector<int>> rot_bins_1rh; // rotation bins at 1rh radius
    vector<vector<int>> rot_bins_3rh; // rotation bins at 3rh radius
    vector<vector<int>> rot_bins_5rh; // rotation bins at 5rh radius
    vector<vector<int>> rot_bins_7rh; // rotation bins at 7rh radius
    vector<vector<int>> rot_bins_rarh; // rotation bins at rarh radius

    vector<shared_ptr<Hetatom>> place() const override {
        // dereference the values we'll need for better performance
        const vector<int> bins = grid->get_bins();
        vector<vector<vector<char>>>& gref = grid->grid;
        const int ra = grid->ra, rh = grid->rh;
        const int r_eff = ra+rh;

        // we define a helper lambda
        vector<shared_ptr<Hetatom>> placed_water;
        auto add_loc = [&] (const vector<int> v) {
            shared_ptr<Hetatom> a = Hetatom::create_new_water(grid->to_xyz(v));
            grid->add(a);
            grid->expand_volume(a);
            placed_water.push_back(a);
        };

        vector<shared_ptr<Atom>> atoms = grid->get_protein_atoms();
        for (auto const& a : atoms) {
            const vector<int>& loc = grid->members.at(a);
            const int x = loc[0], y = loc[1], z = loc[2];

            cout << "(" << x << ", " << y << ", " << z << ")" << endl;


            for (int i = 0; i < rot_bins_rarh.size(); i++) {
                int xr = x + rot_bins_rarh[i][0], yr = y + rot_bins_rarh[i][1], zr = z + rot_bins_rarh[i][2]; // new coordinates
                
                // check bounds
                if (xr < 0) xr = 0;
                if (xr >= bins[0]) xr = bins[0]-1;
                if (yr < 0) yr = 0;
                if (yr >= bins[1]) yr = bins[1]-1;
                if (zr < 0) zr = 0;
                if (zr >= bins[2]) zr = bins[2]-1;

                // we have to make sure we don't check the direction of the atom we are trying to place this water on
                const vector<int> skip_bin = {xr-rot_bins_1rh[i][0], yr-rot_bins_1rh[i][1], zr-rot_bins_1rh[i][2]};
                if (gref[xr][yr][zr] == 0 && collision_check({xr, yr, zr}, skip_bin)) {add_loc({xr, yr, zr});};
            }
        }

        return placed_water;
    }

private:
    /**
     * @brief Check if a water molecule can be placed at the given location. 
     * @param loc the location to be checked. 
     * @param origin location to be excluded from
     * @return True if this is an acceptable location, false otherwise.
     */
    bool collision_check(const vector<int> loc, const vector<int> skip_bin) const {
        // dereference the values we'll need for better performance
        vector<vector<vector<char>>>& gref = grid->grid;
        const vector<int> bins = grid->get_bins();
        const int ra = grid->ra, rh = grid->rh;

        // check for collisions at 1rh
        for (auto const& rot : rot_bins_1rh) {
            int xr = loc[0] + rot[0], yr = loc[1] + rot[1], zr = loc[2] + rot[2]; // new coordinates
            if (vector({xr, yr, zr}) == skip_bin) {continue;}

            if (xr < 0) xr = 0;
            if (xr >= bins[0]) xr = bins[0]-1;
            if (yr < 0) yr = 0;
            if (yr >= bins[1]) yr = bins[1]-1;
            if (zr < 0) zr = 0;
            if (zr >= bins[2]) zr = bins[2]-1;

            if (gref[xr][yr][zr] != 0) {
                // cout << "Collision1 at (" << xr << ", " << yr << ", " << zr << ")" << endl;
                // cout << gref[xr][yr][zr] << endl;
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