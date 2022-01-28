// includes
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <list>

// my own includes
#include "hydrate/Grid.h"
#include "data/Atom.h"
#include "data/Body.h"
#include "data/Hetatom.h"
#include "hydrate/AxesPlacement.cpp"
#include "hydrate/RadialPlacement.cpp"
#include "hydrate/JanPlacement.cpp"
#include "hydrate/CounterCulling.cpp"
#include "hydrate/OutlierCulling.cpp"
#include "settings.h"
#include "math/Vector3.h"

using boost::format;
using std::vector, std::string, std::cout, std::endl, std::shared_ptr, std::unique_ptr;
using namespace setting::grid;

Grid::Grid(const Axis3D axes, double width, double ra, double rh, PlacementStrategyChoice psc, CullingStrategyChoice csc) : axes(axes) {
    long long int total_bins = (long long) axes.x.bins*axes.y.bins*axes.z.bins;
    if (total_bins > 32e9) {
        throw except::invalid_argument("Error in Grid: Too many bins.");
    } else if (total_bins > 4e9) {
        print_err("Warning in Grid: Consider lowering the number of bins.");
    }
    this->width = width;
    this->grid = vector(axes.x.bins, vector<vector<char>>(axes.y.bins, vector<char>(axes.z.bins, 0)));
    this->set_radius_atoms(ra);
    this->set_radius_water(rh);

    switch (psc) {
        case AxesStrategy: 
            water_placer = std::make_unique<AxesPlacement>(this);
            break;
        case RadialStrategy:
            water_placer = std::make_unique<RadialPlacement>(this);
            break;
        case JanStrategy: 
            water_placer = std::make_unique<JanPlacement>(this);
            break;
        default: 
            throw except::unknown_argument("Error in Grid::Grid: Unkown PlacementStrategy");
    }

    switch (csc) {
        case CounterStrategy: 
            water_culler = std::make_unique<CounterCulling>(this);
            break;
        case OutlierStrategy: 
            water_culler = std::make_unique<OutlierCulling>(this);
            break;
        default: 
            throw except::unknown_argument("Error in Grid::Grid: Unkown CullingStrategy");
    }
}

vector<Hetatom> Grid::hydrate() {
    vector<GridMember<Hetatom>> placed_water = find_free_locs(); // the molecules which were placed by the find_free_locs method
    water_culler->set_target_count(setting::grid::percent_water*a_members.size()); // target is 10% of atoms
    auto temp = water_culler->cull(placed_water);
    return temp;
}

vector<GridMember<Hetatom>> Grid::find_free_locs() {
    // a quick check to verify there are no water molecules already present
    if (w_members.size() != 0) {
            print_err("Warning in Grid::find_free_locs: Attempting to hydrate a grid which already contains water!");
    }
    expand_volume();

    // place the water molecules with the chosen strategy
    vector<GridMember<Hetatom>> placed_water = water_placer->place();
    return placed_water;
}

vector<vector<int>> Grid::bounding_box() const {
    if (__builtin_expect(a_members.size() == 0, false)) {
        throw except::invalid_operation("Error in Grid::bounding_box: Calculating a boundary box for a grid with no members!");
    }

    // initialize the bounds as large as possible
    vector<vector<int>> box = {{axes.x.bins, 0}, {axes.y.bins, 0}, {axes.z.bins, 0}};
    for (const auto& atom : a_members) {
        for (int i = 0; i < 3; i++) {
            if (box[i][0] > atom.loc[i]) box[i][0] = atom.loc[i]; // min
            if (box[i][1] < atom.loc[i]) box[i][1] = atom.loc[i]+1; // max. +1 since this will often be used as loop limits
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

vector<Hetatom> Grid::get_waters() const {
    vector<Hetatom> atoms(w_members.size());
    int i = 0; // counter
    for (const auto& water : w_members) {
        atoms[i] = water.atom;
        i++;
    }
    atoms.resize(i);
    return atoms;
}

vector<Atom> Grid::get_atoms() const {
    vector<Atom> atoms(a_members.size());
    int i = 0; // counter
    for (const auto& atom : a_members) {
        atoms[i] = atom.atom;
        i++;
    }
    atoms.resize(i);
    return atoms;
}

void Grid::expand_volume() {
    // iterate through each member location
    for (auto& atom : a_members) {
        expand_volume(atom);
    }

    for (auto& water : w_members) {
        expand_volume(water);
    }
}

void Grid::expand_volume(GridMember<Atom>& atom) {
    if (atom.expanded_volume) {return;} // check if this location has already been expanded

    atom.expanded_volume = true; // mark this location as expanded
    expand_volume(atom.loc, false); // do the expansion
}

void Grid::expand_volume(GridMember<Hetatom>& water) {
    if (water.expanded_volume) {return;} // check if this location has already been expanded

    water.expanded_volume = true; // mark this location as expanded
    expand_volume(water.loc, true); // do the expansion
}

void Grid::expand_volume(const vector<int>& loc, const bool is_water) {
    char marker = is_water ? 'h' : 'a';
    const int x = loc[0], y = loc[1], z = loc[2];

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int r = is_water ? rh : ra; // determine which radius to use for the expansion
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, axes.x.bins); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, axes.y.bins); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, axes.z.bins); // zminus and zplus

    // loop over each bin in the box
    int added_volume = 0;
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(x - i, 2) + std::pow(y - j, 2) + std::pow(z - k, 2)) <= r) {
                    if (grid[i][j][k] != 0) {continue;} // skip if the bin is already occupied
                    added_volume++;
                    grid[i][j][k] = marker;
                }
            }
        }
    }

    if (!is_water) {volume += added_volume;};
}

vector<GridMember<Atom>> Grid::add(const Body* const body) {
    return add(body->protein_atoms);
}

void Grid::remove(const Body* const body) {
    remove(body->protein_atoms);
}

GridMember<Atom> Grid::add(const Atom& atom, const bool expand) {
    vector<int> loc = to_bins(atom.coords);
    const int &x = loc[0], &y = loc[1], &z = loc[2];

    // sanity check
    const bool out_of_bounds = x >= axes.x.bins || y >= axes.y.bins || z >= axes.z.bins || x+y+z < 0;
    if (__builtin_expect(out_of_bounds, false)) {
        throw except::out_of_bounds("Error in Grid::add: Atom is located outside the grid!\nLocation: " + atom.coords.to_string());
    }

    if (grid[x][y][z] == 0) {volume++;} // can probably be removed

    GridMember gm(atom, loc);
    if (expand) {expand_volume(gm);}
    a_members.push_back(gm);
    grid[x][y][z] = 'A';

    return gm;
}

GridMember<Hetatom> Grid::add(const Hetatom& water, const bool expand) {
    vector<int> loc = to_bins(water.coords);
    const int &x = loc[0], &y = loc[1], &z = loc[2];

    // sanity check
    const bool out_of_bounds = x >= axes.x.bins || y >= axes.y.bins || z >= axes.z.bins || x+y+z < 0;
    if (__builtin_expect(out_of_bounds, false)) {
        throw except::out_of_bounds("Error in Grid::add: Atom is located outside the grid!\nLocation: " + water.coords.to_string());
    }

    GridMember gm(water, loc);
    if (expand) {expand_volume(gm);}
    w_members.push_back(gm);
    grid[x][y][z] = 'H';

    return gm;
}

void Grid::remove(const Atom& atom) {
    auto pos = std::find(a_members.begin(), a_members.end(), atom);
    if (__builtin_expect(pos == a_members.end(), false)) {
        throw except::invalid_operation("Error in Grid::remove: Attempting to remove an atom which is not part of the grid!");
    }

    GridMember<Atom>& member = *pos;
    const int x = member.loc[0], y = member.loc[1], z = member.loc[2];

    deflate_volume(member);
    a_members.erase(pos);
    grid[x][y][z] = 0;
    volume--;
}

void Grid::remove(const Hetatom& water) {
    auto pos = std::find(w_members.begin(), w_members.end(), water);
    if (__builtin_expect(pos == w_members.end(), false)) {
        throw except::invalid_operation("Error in Grid::remove: Attempting to remove an atom which is not part of the grid!");
    }

    GridMember<Hetatom>& member = *pos;
    const int x = member.loc[0], y = member.loc[1], z = member.loc[2];

    deflate_volume(member);
    w_members.erase(pos);
    grid[x][y][z] = 0;
}

void Grid::remove(const vector<Atom>& atoms) {
    std::cout << "REMOVE 1" << std::endl;    
    
    // we make a vector of all possible uids
    // std::unordered_map<int, int> removed;
    vector<bool> removed(Atom::uid_counter);
    // and fill it with the uids that should be removed
    std::for_each(atoms.begin(), atoms.end(), [&removed] (const Atom& water) {removed[water.uid] = true;});

    std::cout << "REMOVE 2" << std::endl;    
    size_t index = 0; // current index in removed_waters
    vector<GridMember<Atom>> removed_atoms(atoms.size()); // the waters which will be removed
    auto predicate = [&removed, &removed_atoms, &index] (const GridMember<Atom>& gm) {
        if (removed[gm.atom.uid]) { // now we can simply look up in our removed vector to determine if an element should be removed
            removed_atoms[index++] = gm;
            return true;
        }
        return false;
    };

    std::cout << "REMOVE 3" << std::endl;    
    // we save the sizes so we can make a sanity check after the removal    
    size_t prev_size = a_members.size();
    a_members.remove_if(predicate);
    size_t cur_size = a_members.size();

    // sanity check
    if (__builtin_expect(prev_size - cur_size != atoms.size(), false)) {
        throw except::invalid_operation("Error in Grid::remove: Something went wrong.");
    }

    std::cout << "REMOVE 4" << std::endl;    
    // clean up the grid
    for (auto& atom : removed_atoms) {
        const int x = atom.loc[0], y = atom.loc[1], z = atom.loc[2];
        deflate_volume(atom);
        grid[x][y][z] = 0;
        volume--;
    }
}

void Grid::remove(const vector<Hetatom>& waters) {
    // we make a vector of all possible uids
    vector<bool> removed(Atom::uid_counter);
    // and fill it with the uids that should be removed
    std::for_each(waters.begin(), waters.end(), [&removed] (const Hetatom& water) {removed[water.uid] = true;});

    size_t index = 0; // current index in removed_waters
    vector<GridMember<Hetatom>> removed_waters(waters.size()); // the waters which will be removed
    auto predicate = [&removed, &removed_waters, &index] (const GridMember<Hetatom>& gm) {
        if (removed[gm.atom.uid]) { // now we can simply look up in our removed vector to determine if an element should be removed
            removed_waters[index++] = gm;
            return true;
        }
        return false;
    };

    // we save the sizes so we can make a sanity check after the removal    
    size_t prev_size = w_members.size();
    w_members.remove_if(predicate);
    size_t cur_size = w_members.size();

    // sanity check
    if (__builtin_expect(prev_size - cur_size != waters.size(), false)) {
        throw except::invalid_operation("Error in Grid::remove: Something went wrong.");
    }

    // clean up the grid
    for (auto& atom : removed_waters) {
        const int x = atom.loc[0], y = atom.loc[1], z = atom.loc[2];
        deflate_volume(atom);
        grid[x][y][z] = 0;
    }
}

void Grid::deflate_volume() {
    // iterate through each member location
    for (auto& atom : a_members) {
        deflate_volume(atom);
    }

    for (auto& water : w_members) {
        deflate_volume(water);
    }
}

void Grid::deflate_volume(GridMember<Atom>& atom) {
    if (!atom.expanded_volume) {return;} // check if this location has already been expanded

    atom.expanded_volume = false; // mark this location as expanded
    deflate_volume(atom.loc, false); // do the expansion
}

void Grid::deflate_volume(GridMember<Hetatom>& water) {
    if (!water.expanded_volume) {return;} // check if this location has already been expanded

    water.expanded_volume = false; // mark this location as expanded
    deflate_volume(water.loc, true); // do the expansion
}

void Grid::deflate_volume(const vector<int>& loc, const bool is_water) {
    const int x = loc[0], y = loc[1], z = loc[2];

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int r = is_water ? rh : ra; // determine which radius to use for the expansion
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, axes.x.bins); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, axes.y.bins); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, axes.z.bins); // zminus and zplus

    // loop over each bin in the box
    int removed_volume = -1; // -1 because we overwrite the center
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(loc[0] - i, 2) + std::pow(loc[1] - j, 2) + std::pow(loc[2] - k, 2)) <= r) {
                    if (grid[i][j][k] == 0) {continue;} // skip if bin is already empty
                    removed_volume++;
                    grid[i][j][k] = 0;
                }
            }
        }
    }

    if (is_water) {
        grid[x][y][z] = 'H'; // replace the center
    } else {
        grid[x][y][z] = 'A';
        volume -= removed_volume; // only the actual atoms contributes to the volume
    }
}

vector<int> Grid::get_bins() const {
    return {axes.x.bins, axes.y.bins, axes.z.bins};
}

vector<int> Grid::to_bins(const Vector3& v) const {
    int binx = std::round((v[0] - axes.x.min)/width);
    int biny = std::round((v[1] - axes.y.min)/width);
    int binz = std::round((v[2] - axes.z.min)/width);
    return {binx, biny, binz};
}

Vector3 Grid::to_xyz(const vector<int>& v) const {
    double x = axes.x.min + width*v[0];
    double y = axes.y.min + width*v[1];
    double z = axes.z.min + width*v[2];
    return {x, y, z};
}

double Grid::get_volume() {
    expand_volume();
    return pow(width, 3)*volume;
}