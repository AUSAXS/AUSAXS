#pragma once

// includes
#include <vector>
#include <map>
#include "boost/format.hpp"
#include <utility>

// ROOT
#include <TVector3.h>

// my own includes
#include "Grid.h"
#include "data/Atom.cpp"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

Grid::Grid(TVector3 base, double width, vector<int> bins, double ra, double rh) {
    long long int total_bins = (long long) bins[0]*bins[1]*bins[2];
    if (total_bins > 32e9) {
        print_err("Error in Grid: Too many bins.");
        exit(1);
    } else if (total_bins > 4e9) {
        print_err("Warning in Grid: Consider lowering the number of bins.");
    }
    this->base = base;
    this->width = width;
    this->bins = bins;
    this->grid = vector(bins[0], vector<vector<char>>(bins[1], vector<char>(bins[2], 0)));
    this->set_radius_atoms(ra);
    this->set_radius_water(rh);
}

void Grid::add(vector<shared_ptr<Atom>>* atoms) {
    for (auto const& a : *atoms) {
        add(a);
    }
}

void Grid::expand_volume() {
    vol_expanded = true;

    // iterate through each member location
    for (const auto& pair : members) {
        expand_volume(pair.first);
    }
}

vector<shared_ptr<Atom>> Grid::hydrate(int reduce = 3) {
    vector<shared_ptr<Atom>> hydration_atoms;
    vector<vector<int>> hydration_slots = find_free_locs();

    int c = 0; // counter
    for (vector<int> v : hydration_slots) {
        c++;
        shared_ptr<Atom> a = Atom::create_new_water(to_xyz(v));
        if (reduce != 0) {
            if (c % reduce != 0) {
                continue;
            }
        }
        add(a);
        expand_volume(*a);
        hydration_atoms.push_back(a);
    }
    return hydration_atoms;
}

vector<vector<int>> Grid::find_free_locs() {
    // a quick check to verify there are no water molecules already present
    for (const auto& pair : members) {
        if (pair.first.is_water()) {
            print_err("Warning in Grid::find_free_locs: Attempting to hydrate a grid which already contains water!");
        }
    }

    if (!vol_expanded) {
        expand_volume();
    }

    vector<vector<int>> bounds = bounding_box();
    // add ra+rh extra space in each direction to account for the volume of both atoms and water molecules
    for (int i = 0; i < 3; i++) {
        bounds[i][0] = std::max(bounds[i][0] - (ra+rh), 0);
        bounds[i][1] = std::min(bounds[i][1] + (ra+rh) + 1, bins[i]); // +1 since this range is inclusive, but the following loop is not
    }

    // we define two helper functions so I can make the checks in the inner loop one-liners
    vector<vector<int>> available_locs;        
    auto check_loc = [&] (const vector<int> v) {return check_collisions(v, available_locs);};
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
    cout << "Found " << available_locs.size() << " available HOH spots." << endl;
    return available_locs;
}

bool Grid::check_collisions(const vector<int> loc, const vector<vector<int>> other_molecules) const {
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

vector<vector<int>> Grid::bounding_box() const {
    if (members.size() == 0) {
        print_err("Error in Grid::bounding_box: Calculating a boundary box for a grid with no members!");
        exit(1);
    }

    // initialize the bounds as large as possible
    vector<vector<int>> box = {{bins[0], 0}, {bins[1], 0}, {bins[2], 0}};
    for (const auto& pair : members) {
        const vector<int>& loc = pair.second;
        for (int i = 0; i < 3; i++) {
            if (box[i][0] > loc[i]) box[i][0] = loc[i]; // min
            if (box[i][1] < loc[i]) box[i][1] = loc[i]; // max
        }
    }
    return box;
}

vector<vector<vector<char>>>* Grid::get_grid() {return &grid;}

void Grid::set_radius_atoms(double radius) {
    int new_r = int(radius/width); // convert the radius to a "bin-radius"
    if (this->ra != 0 && this->ra != new_r) {
        print_err("Warning in Grid::set_radius: The radius is already set for this grid!");
    }
    this->ra = new_r;
}

void Grid::set_radius_water(double radius) {
    int new_r = int(radius/width); // convert the radius to a "bin-radius"
    if (this->rh != 0 && this->rh != new_r) {
        print_err("Warning in Grid::set_radius: The radius is already set for this grid!");
    }
    this->rh = new_r;
}

vector<Atom*> Grid::get_hydration_atoms() const {
    vector<Atom*> atoms;
    for (const auto& pair : members) {
        Atom a = pair.first;
        if (a.is_water()) {
            atoms.push_back(&a);
        }
    }
    return atoms;
}

double Grid::get_volume() {return pow(width, 3)*volume;}

void Grid::expand_volume(const Atom atom) {
    vector<int> loc = members.at(atom);
    int r = atom.is_water() ? rh : ra; // determine which radius to use for the expansion
    char marker = atom.is_water() ? 'H' : 'A';

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    vector<vector<int>> bounds(3, vector<int>(2, 0));
    for (int i = 0; i < 3; i++) {
        bounds[i][0] = std::max(loc[i] - r, 0);
        bounds[i][1] = std::min(loc[i] + r + 1, bins[i]); // +1 since this range is inclusive, while the following for-loop is not
    }

    // loop over each bin in the box
    for (int i = bounds[0][0]; i < bounds[0][1]; i++) {
        for (int j = bounds[1][0]; j < bounds[1][1]; j++) {
            for (int k = bounds[2][0]; k < bounds[2][1]; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(loc[0] - i, 2) + std::pow(loc[1] - j, 2) + std::pow(loc[2] - k, 2)) <= r) {
                    volume++;
                    grid[i][j][k] = marker;
                }
            }
        }
    }
}

void Grid::add(shared_ptr<Atom> atom) {
    volume++;
    vector<int> loc = to_bins(atom->get_coords());

    // sanity check
    if (loc[0] >= bins[0] || loc[1] >= bins[1] || loc[2] >= bins[2]) {
        print_err("Error in Grid::add: Atom is located outside the grid!");
        exit(1);
    }
    members.insert({*atom, loc});
    grid[loc[0]][loc[1]][loc[2]] = atom->is_water() ? 'H' : 'A';
}

vector<int> Grid::to_bins(TVector3 v) {
    int binx = std::round((v[0] - base.X())/width);
    int biny = std::round((v[1] - base.Y())/width);
    int binz = std::round((v[2] - base.Z())/width);
    return {binx, biny, binz};
}

TVector3 Grid::to_xyz(vector<int> v) {
    double x = base.X() + width*v[0];
    double y = base.Y() + width*v[1];
    double z = base.Z() + width*v[2];
    return {x, y, z};
}
