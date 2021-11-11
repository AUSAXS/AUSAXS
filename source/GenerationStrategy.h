#pragma once

// includes
#include <TVector3.h>
#include "data/Atom.cpp"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

class PlacementStrategy {
public:
    PlacementStrategy(vector<vector<vector<char>>>& grid, vector<int> bins, int ra, int rh){
        this->grid = grid;
        this->ra = ra;
        this->rh = rh;
    }

    virtual void place_water(const vector<int> loc) = 0;

protected: 
    vector<vector<vector<char>>> grid;
    vector<int> bins;
    int ra, rh;

private:
    virtual bool collision_check(const vector<int> loc) const = 0;

};

class ClassicPlacement : public PlacementStrategy {
public:
    ClassicPlacement(vector<vector<vector<char>>>& grid, vector<int> bins, int ra, int rh) : PlacementStrategy(grid, bins, ra, rh) {}

    void place_water(vector<vector<int>> bounds) {
        // we define two helper functions so I can make the checks in the inner loop one-liners
        vector<vector<int>> available_locs;
        auto check_loc = [&] (const vector<int> v) {return collision_check(v);};
        auto add_loc = [&] (const vector<int> v) {
            grid[v[0]][v[1]][v[2]] = 'H';
            available_locs.push_back({v[0], v[1], v[2]});
        };

        // loop over the minimum bounding box as found above
        for (int i = bounds[0][0]; i < bounds[0][1]; i++) {
            for (int j = bounds[1][0]; j < bounds[1][1]; j++) {
                for (int k = bounds[2][0]; k < bounds[2][1]; k++) {
                    // if this spot is part of an atom
                    if (grid[i][j][k] == 'A') {
                        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
                        int xmin = std::max(i-rh, 0), xmax = std::min(i+rh, bins[0]);
                        int ymin = std::max(j-rh, 0), ymax = std::min(j+rh, bins[1]);
                        int zmin = std::max(k-rh, 0), zmax = std::min(k+rh, bins[2]);

                        // check collisions for x ± r_eff
                        if ((grid[xmin][j][k] == 0) && check_loc({xmin, j, k})) add_loc({xmin, j, k});
                        if ((grid[xmax][j][k] == 0) && check_loc({xmax, j, k})) add_loc({xmax, j, k});

                        // check collisions for y ± r_eff
                        if ((grid[i][ymin][k] == 0) && check_loc({i, ymin, k})) add_loc({i, ymin, k});
                        if ((grid[i][ymax][k] == 0) && check_loc({i, ymax, k})) add_loc({i, ymax, k});

                        // check collisions for z ± r_eff
                        if ((grid[i][j][zmin] == 0) && check_loc({i, j, zmin})) add_loc({i, j, zmin});
                        if ((grid[i][j][zmax] == 0) && check_loc({i, j, zmax})) add_loc({i, j, zmax});
                    }
                }
            }
        }
    }

private:
    bool collision_check(const vector<int> loc) const {
        const int x = loc[0], y = loc[1], z = loc[2];

        // check collision with other atoms
        int r_eff = ra + rh; // effective radius for water/atom collisions
        return true;
    }
};