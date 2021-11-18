// includes
#include <vector>
#include <map>
#include "boost/format.hpp"
#include <utility>

// ROOT
#include <TVector3.h>

// my own includes
#include "Grid.h"
#include "data/Atom.h"
#include "data/Hetatom.h"
#include "AxesPlacement.cpp"
#include "RadialPlacement.cpp"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

Grid::Grid(TVector3 base, double width, vector<int> bins, double ra, double rh, PlacementStrategyChoice ch = AxesStrategy) {
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

    if (ch == AxesStrategy) {water_placer = std::make_unique<AxesPlacement>(this);}
    else if (ch == RadialStrategy) {water_placer = std::make_unique<RadialPlacement>(this);}
    else {print_err("Error in Grid::Grid: Unkown PlacementStrategy");}
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

vector<shared_ptr<Hetatom>> Grid::hydrate(int reduce = 3) {
    vector<shared_ptr<Hetatom>> placed_water = find_free_locs();

    int c = 0; // counter
    for (const auto& a : placed_water) {
        if (reduce != 0) {
            if (c % reduce != 0) {
                remove(a);
                placed_water.erase(placed_water.begin()+c);
                continue;
            }
        }
        c++;
    }
    return placed_water;
}

vector<shared_ptr<Hetatom>> Grid::find_free_locs() {
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

    // place the water molecules with the chosen strategy
    vector<shared_ptr<Hetatom>> placed_water = water_placer->place(bounds);

    cout << "Placed " << placed_water.size() << " HOH molecules." << endl;
    return placed_water;
}

// bool Grid::check_collisions(const vector<int> loc, const vector<vector<int>> other_molecules) const {
//     const int x = loc[0], y = loc[1], z = loc[2];

//     // check collision with other atoms
//     int r_eff = ra + rh; // effective radius for water/atom collisions
//     for (auto const& pair : members) {
//         const vector<int>& v = pair.second;

//         // check if the point is inside the spherical volume of loc
//         if (sqrt(pow(x-v[0], 2) + pow(y-v[1], 2) + pow(z-v[2], 2)) < r_eff) {
//             return false;
//         }
//     }

//     // check collision with other water molecules
//     r_eff = 2*rh; // effective radius for water/water collisions
//     for (auto const& a : other_molecules) {
//         // check if the point is inside the spherical volume of loc
//         if (sqrt(pow(x-a[0], 2) + pow(y-a[1], 2) + pow(z-a[2], 2)) < r_eff) {
//             return false;
//         }
//     }
//     return true;
// }

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

void Grid::expand_volume(const Atom atom) {
    vector<int> loc = members.at(atom);
    char marker = atom.is_water() ? 'h' : 'a';

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int r = atom.is_water() ? rh : ra; // determine which radius to use for the expansion
    int xm = std::max(loc[0]-r, 0), xp = std::min(loc[0]+r+1, bins[0]); // xminus and xplus
    int ym = std::max(loc[1]-r, 0), yp = std::min(loc[1]+r+1, bins[1]); // yminus and yplus
    int zm = std::max(loc[2]-r, 0), zp = std::min(loc[2]+r+1, bins[2]); // zminus and zplus

    // loop over each bin in the box
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(loc[0] - i, 2) + std::pow(loc[1] - j, 2) + std::pow(loc[2] - k, 2)) <= r) {
                    volume++;
                    grid[i][j][k] = marker;
                }
            }
        }
    }

    grid[loc[0]][loc[1]][loc[2]] = atom.is_water() ? 'H' : 'A'; // replace the center with a capital letter (better than doing another if-statement in the loop)
}

void Grid::remove(shared_ptr<Atom> atom) {
    if (members.count(*atom) == 0) {
        print_err("Error in Grid::remove: Attempting to remove an atom which is not part of the grid!");
        exit(1);
    }

    vector<int> loc = members.at(*atom);
    members.erase(*atom);
    grid[loc[0]][loc[1]][loc[2]] = 0;

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int r = atom->is_water() ? rh : ra; // determine which radius to use for the expansion
    int xm = std::max(loc[0]-r, 0), xp = std::min(loc[0]+r+1, bins[0]); // xminus and xplus
    int ym = std::max(loc[1]-r, 0), yp = std::min(loc[1]+r+1, bins[1]); // yminus and yplus
    int zm = std::max(loc[2]-r, 0), zp = std::min(loc[2]+r+1, bins[2]); // zminus and zplus

    // loop over each bin in the box
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(loc[0] - i, 2) + std::pow(loc[1] - j, 2) + std::pow(loc[2] - k, 2)) <= r) {
                    volume--;
                    grid[i][j][k] = 0;
                }
            }
        }
    }
}

void Grid::add(shared_ptr<Atom> atom) {
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
