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
#include "CounterCulling.cpp"
#include "OutlierCulling.cpp"
#include "settings.h"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::shared_ptr, std::unique_ptr;
using namespace ROOT;
using namespace setting::grid;

Grid::Grid(TVector3 base, double width, vector<int> bins, double ra, double rh, PlacementStrategyChoice psc, CullingStrategyChoice csc) {
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

    switch (psc) {
        case AxesStrategy: 
            water_placer = std::make_unique<AxesPlacement>(this);
            break;
        case RadialStrategy:
            water_placer = std::make_unique<RadialPlacement>(this);
            break;
        default: 
            print_err("Error in Grid::Grid: Unkown PlacementStrategy");
            exit(1);
    }

    switch (csc) {
        case CounterStrategy: 
            water_culler = std::make_unique<CounterCulling>(this);
            break;
        case OutlierStrategy: 
            water_culler = std::make_unique<OutlierCulling>(this);
            break;
        default: 
            print_err("Error in Grid::Grid: Unkown CullingStrategy");
            exit(1);        
    }
}

void Grid::expand_volume() {
    vol_expanded = true;

    // iterate through each member location
    for (const auto& pair : members) {
        expand_volume(pair.first);
    }
}

vector<shared_ptr<Hetatom>> Grid::hydrate() {
    vector<shared_ptr<Hetatom>> placed_water = find_free_locs(); // the molecules which were placed by the find_free_locs method
    water_culler->set_target_count(setting::grid::percent_water*(members.size()-placed_water.size())); // target is 10% of atoms
    return water_culler->cull(placed_water);
}

vector<shared_ptr<Hetatom>> Grid::find_free_locs() {
    // a quick check to verify there are no water molecules already present
    for (const auto& pair : members) {
        if (pair.first->is_water()) {
            print_err("Warning in Grid::find_free_locs: Attempting to hydrate a grid which already contains water!");
        }
    }

    if (!vol_expanded) {
        expand_volume();
    }

    // place the water molecules with the chosen strategy
    vector<shared_ptr<Hetatom>> placed_water = water_placer->place();
    cout << "Placed " << placed_water.size() << " HOH molecules." << endl;
    return placed_water;
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

vector<shared_ptr<Atom>> Grid::get_hydration_atoms() const {
    vector<shared_ptr<Atom>> atoms(members.size());
    int i = 0; // counter
    for (const auto& pair : members) {
        if (pair.first->is_water()) {
            atoms[i] = pair.first;
            i++;
        }
    }
    atoms.resize(i);
    return atoms;
}

vector<shared_ptr<Atom>> Grid::get_protein_atoms() const {
    vector<shared_ptr<Atom>> atoms(members.size());
    int i = 0; // counter
    for (const auto& pair : members) {
        if (!pair.first->is_water()) {
            atoms[i] = pair.first;
            i++;
        }
    }
    atoms.resize(i);
    return atoms;
}

void Grid::expand_volume(const shared_ptr<Atom> atom) {
    vector<int> loc = members.at(atom);
    char marker = atom->is_water() ? 'h' : 'a';

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int r = atom->is_water() ? rh : ra; // determine which radius to use for the expansion
    int xm = std::max(loc[0]-r, 0), xp = std::min(loc[0]+r+1, bins[0]-1); // xminus and xplus
    int ym = std::max(loc[1]-r, 0), yp = std::min(loc[1]+r+1, bins[1]-1); // yminus and yplus
    int zm = std::max(loc[2]-r, 0), zp = std::min(loc[2]+r+1, bins[2]-1); // zminus and zplus

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

    grid[loc[0]][loc[1]][loc[2]] = atom->is_water() ? 'H' : 'A'; // replace the center with a capital letter (better than doing another if-statement in the loop)
}

void Grid::remove(const shared_ptr<Atom> atom) {
    if (members.count(atom) == 0) {
        print_err("Error in Grid::remove: Attempting to remove an atom which is not part of the grid!");
        exit(1);
    }

    vector<int> loc = members.at(atom);
    const int x = loc[0], y = loc[1], z = loc[2];
    members.erase(atom);
    grid[x][y][z] = 0;

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int r = atom->is_water() ? rh : ra; // determine which radius to use for the expansion
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, bins[0]-1); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, bins[1]-1); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, bins[2]-1); // zminus and zplus

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

void Grid::add(const shared_ptr<Atom> atom) {
    vector<int> loc = to_bins(atom->get_coords());
    const int x = loc[0], y = loc[1], z = loc[2];

    // sanity check
    if (x >= bins[0] || y >= bins[1] || z >= bins[2]) {
        print_err("Error in Grid::add: Atom is located outside the grid!");
        TVector3 coords = atom->get_coords();
        print_err((format("Location: (%1%, %2%, %3%)") % coords[0] % coords[1] % coords[2]).str());
        exit(1);
    }
    members.insert({atom, loc});
    grid[x][y][z] = atom->is_water() ? 'H' : 'A';
}

vector<int> Grid::to_bins(const TVector3& v) const {
    int binx = std::round((v[0] - base.X())/width);
    int biny = std::round((v[1] - base.Y())/width);
    int binz = std::round((v[2] - base.Z())/width);
    return {binx, biny, binz};
}

TVector3 Grid::to_xyz(const vector<int>& v) const {
    double x = base.X() + width*v[0];
    double y = base.Y() + width*v[1];
    double z = base.Z() + width*v[2];
    return {x, y, z};
}
