#pragma once

#include "Grid.cpp"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

class Grid_classic : public Grid {
private:
    /**
     * @brief Check if a water molecule can be placed at the given location. 
     *        This checks collisions with both other water molecules and other atoms. 
     * @param loc the location to be checked. 
     * @param other_molecules the other water molecules which have already been placed.
     * @return True if this is an acceptable location, false otherwise.
     */
    bool check_collisions(const vector<int> loc, const vector<vector<int>> other_molecules) const override {
        const int x = loc[0], y = loc[1], z = loc[2];

        double amin = 100, wmin = 100;

        // check collision with other atoms
        int r_eff = ra + rh; // effective radius for water/atom collisions
        for (auto const& pair : members) {
            const vector<int>& v = pair.second;

            // check if the point is inside the spherical volume of loc
            double d = sqrt(pow(x-v[0], 2) + pow(y-v[1], 2) + pow(z-v[2], 2));
            amin = std::min(amin, d);
            if (d < r_eff) {
                // cout << format("Atom: (%1%, %2%, %3%), water: (%4%, %5%, %6%)") % v[0] % v[1] % v[2] % x % y % z << endl;
                // cout << "Separation: " << sqrt(pow(x-v[0], 2) + pow(y-v[1], 2) + pow(z-v[2], 2)) << endl << endl;
                return false;
            }
        }

        // check collision with other water molecules
        r_eff = 2*rh; // effective radius for water/water collisions
        for (auto const& a : other_molecules) {
            // check if the point is inside the spherical volume of loc
            double d = sqrt(pow(x-a[0], 2) + pow(y-a[1], 2) + pow(z-a[2], 2));
            wmin = std::min(wmin, d);
            if (sqrt(pow(x-a[0], 2) + pow(y-a[1], 2) + pow(z-a[2], 2)) < r_eff) {
                return false;
            }            
        }
        // cout << "Water molecule placed! (amin, wmin): (" << amin << ", " << wmin << ")" << endl;
        return true;
    }
};