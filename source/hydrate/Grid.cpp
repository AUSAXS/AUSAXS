#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <list>
#include <limits>
#include <typeinfo>
#include <cassert>

#include <hydrate/Grid.h>
#include <data/Atom.h>
#include <data/Body.h>
#include <data/Hetatom.h>
#include <hydrate/AxesPlacement.h>
#include <hydrate/RadialPlacement.h>
#include <hydrate/JanPlacement.h>
#include <hydrate/CounterCulling.h>
#include <hydrate/OutlierCulling.h>
#include <settings.h>
#include <math/Vector3.h>

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace setting::grid;
using namespace grid;

Grid::Grid(const Axis3D& axes, double width, double ra, double rh, PlacementStrategyChoice psc, CullingStrategyChoice csc) : axes(axes) {
    setup(width, ra, rh, psc, csc);
}

Grid::Grid(const vector<Atom>& atoms, double width, double ra, double rh, PlacementStrategyChoice psc, CullingStrategyChoice csc) {
    // find the bounding box
    vector<int> imin, imax;
    auto[min, max] = bounding_box(atoms);

    // expand the box by 10%
    for (auto& v : min) {
        if (v < 0) {imin.push_back(std::round(v*(1 + setting::grid::scaling)));}     // if v is smaller than 0, multiply by 1+s
        else {imin.push_back(std::round(v*(1 - setting::grid::scaling)));}           //                    else multiply by 1-s
    }
    for (auto& v : max) {
        if (v > 0) {imax.push_back(std::round(v*(1 + setting::grid::scaling)) + 1);} // if v is larger than 0, multiply by 1+s
        else {imax.push_back(std::round(v*(1 - setting::grid::scaling)) + 1);}       //                   else multiply by 1-s
    }

    // setup the rest of the class members
    axes = Axis3D(imin, imax, setting::grid::width);
    setup(width, ra, rh, psc, csc);

    // finally add the atoms to the grid
    add(atoms);
}

Grid::Grid(const vector<Body>& bodies, double width, double ra, double rh, setting::grid::PlacementStrategyChoice psc, setting::grid::CullingStrategyChoice csc) {
    vector<int> imin, imax;

    // find the total bounding box containing all bodies
    Vector3 min{0, 0, 0}, max{0, 0, 0};
    for (const Body& body : bodies) {
        auto[cmin, cmax] = bounding_box(body.protein_atoms);

        for (int i = 0; i < 3; i++) {
            if (cmin[i] < min[i]) {min[i] = cmin[i];}
            if (cmax[i] > max[i]) {max[i] = cmax[i];}
        }
    }

    // expand the box by 10%
    for (auto& v : min) {
        if (v < 0) {imin.push_back(std::round(v*(1 + setting::grid::scaling)));}
        else {imin.push_back(std::round(v*(1 - setting::grid::scaling)));}
    }
    for (auto& v : max) {
        if (v > 0) {imax.push_back(std::round(v*(1 + setting::grid::scaling)) + 1);}
        else {imax.push_back(std::round(v*(1 - setting::grid::scaling)) + 1);}
    }

    // setup the rest of the class members
    axes = Axis3D(imin, imax, setting::grid::width);
    setup(width, ra, rh, psc, csc);

    // finally add all atoms to the grid
    for (const Body& body : bodies) {
        add(body.protein_atoms);
    }
}

Grid::Grid(const Grid& grid) : Grid(grid.axes, grid.width, grid.ra, grid.rh, setting::grid::psc, setting::grid::csc) {
    this->grid = grid.grid;
    this->volume = grid.volume;
    this->ra = grid.ra;
    this->rh = grid.rh;
    this->a_members = grid.a_members;
    this->w_members = grid.w_members;
}

void Grid::setup(double width, double ra, double rh, PlacementStrategyChoice psc, CullingStrategyChoice csc) {
    long long int total_bins = (long long) axes.x.bins*axes.y.bins*axes.z.bins;
    if (total_bins > 32e9) {
        throw except::invalid_argument("Error in Grid: Too many bins.");
    } else if (total_bins > 4e9) {
        print_err("Warning in Grid: Consider lowering the number of bins.");
    }
    this->width = width;
    this->grid = vector(axes.x.bins, vector(axes.y.bins, vector<char>(axes.z.bins, 0)));
    this->set_radius_atoms(ra);
    this->set_radius_water(rh);

    switch (psc) {
        case PlacementStrategyChoice::AxesStrategy: 
            water_placer = std::make_unique<AxesPlacement>(this);
            break;
        case PlacementStrategyChoice::RadialStrategy:
            water_placer = std::make_unique<RadialPlacement>(this);
            break;
        case PlacementStrategyChoice::JanStrategy: 
            water_placer = std::make_unique<JanPlacement>(this);
            break;
        default: 
            throw except::unknown_argument("Error in Grid::Grid: Unkown PlacementStrategy");
    }

    switch (csc) {
        case CullingStrategyChoice::CounterStrategy: 
            water_culler = std::make_unique<CounterCulling>(this);
            break;
        case CullingStrategyChoice::OutlierStrategy: 
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
    vector<vector<int>> box = {{int(axes.x.bins), 0}, {int(axes.y.bins), 0}, {int(axes.z.bins), 0}};
    for (const auto& atom : a_members) {
        for (int i = 0; i < 3; i++) {
            if (box[i][0] > atom.loc[i]) box[i][0] = atom.loc[i]; // min
            if (box[i][1] < atom.loc[i]) box[i][1] = atom.loc[i]+1; // max. +1 since this will often be used as loop limits
        }
    }
    return box;
}

std::pair<Vector3, Vector3> Grid::bounding_box(const vector<Atom>& atoms) {
    if (__builtin_expect(atoms.size() == 0, false)) {
        throw except::invalid_operation("Error in Grid::bounding_box: Calculating a boundary box for a grid with no members!");
    }

    // initialize the bounds as large as possible
    Vector3 min = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vector3 max = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    for (const auto& atom : atoms) {
        for (int i = 0; i < 3; i++) {
            min[i] = std::min(min[i], atom.coords[i]);
            max[i] = std::max(max[i], atom.coords[i]);
        }
    }
    return std::make_pair(min, max);
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
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, int(axes.x.bins)); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, int(axes.y.bins)); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, int(axes.z.bins)); // zminus and zplus

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
    const bool out_of_bounds = x >= int(axes.x.bins) || y >= int(axes.y.bins) || z >= int(axes.z.bins) || x+y+z < 0;
    if (__builtin_expect(out_of_bounds, false)) {
        throw except::out_of_bounds("Error in Grid::add: Atom is located outside the grid!\nLocation: " + atom.coords.to_string() + "\nBounds: " + axes.to_string());
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
    const bool out_of_bounds = x >= int(axes.x.bins) || y >= int(axes.y.bins) || z >= int(axes.z.bins) || x+y+z < 0;
    if (__builtin_expect(out_of_bounds, false)) {
        throw except::out_of_bounds("Error in Grid::add: Atom is located outside the grid!\nLocation: " + water.coords.to_string() + "\nBounds: " + axes.to_string());
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
    // we make a vector of all possible uids
    std::unordered_map<int, bool> removed;
    // and fill it with the uids that should be removed
    std::for_each(atoms.begin(), atoms.end(), [&removed] (const Atom& atom) {removed[atom.uid] = true;});

    size_t index = 0; // current index in removed_atoms
    vector<GridMember<Atom>> removed_atoms(atoms.size()); // the atoms which will be removed
    auto predicate = [&removed, &removed_atoms, &index] (const GridMember<Atom>& gm) {
        if (removed[gm.atom.uid]) { // now we can simply look up in our removed vector to determine if an element should be removed
            removed_atoms[index++] = gm;
            return true;
        }
        return false;
    };

    // we save the sizes so we can make a sanity check after the removal    
    size_t prev_size = a_members.size();
    a_members.remove_if(predicate);
    size_t cur_size = a_members.size();

    // sanity check
    if (__builtin_expect(prev_size - cur_size != atoms.size(), false)) {
        throw except::invalid_operation("Error in Grid::remove: Expected to remove " + std::to_string(atoms.size()) + " elements, but only " + std::to_string(prev_size - cur_size) + " were actually removed.");
    }

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
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, int(axes.x.bins)); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, int(axes.y.bins)); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, int(axes.z.bins)); // zminus and zplus

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

void Grid::clear_waters() {
    vector<Hetatom> waters;
    waters.reserve(w_members.size());

    std::for_each(w_members.begin(), w_members.end(), [&waters] (const GridMember<Hetatom>& water) {waters.push_back(water.atom);});
    remove(waters);
    assert(w_members.size() == 0);
}

vector<int> Grid::get_bins() const {
    return {int(axes.x.bins), int(axes.y.bins), int(axes.z.bins)};
}

vector<int> Grid::to_bins(const Vector3& v) const {
    int binx = std::round((v.x() - axes.x.min)/width);
    int biny = std::round((v.y() - axes.y.min)/width);
    int binz = std::round((v.z() - axes.z.min)/width);
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

Grid Grid::copy() const {
    return Grid(*this);
}

Grid& Grid::operator=(const Grid& rhs) {
    grid = rhs.grid;
    a_members = rhs.a_members;
    w_members = rhs.w_members;
    width = width;
    volume = rhs.volume;
    ra = rhs.ra;
    rh = rhs.rh;
    axes = rhs.axes;
    // culler & placer cannot be modified after program is run, so they'll automatically be equal always
    return *this;
}

bool Grid::operator==(const Grid& rhs) const {
    // we do everything but check the contents of the grid. 
    if (volume != rhs.volume) {return false;}
    if (a_members.size() != rhs.a_members.size()) {return false;}
    if (w_members.size() != rhs.w_members.size()) {return false;}
    if (ra != rhs.ra) {return false;}
    if (rh != rhs.rh) {return false;}
    if (typeid(water_culler) != typeid(rhs.water_culler)) {return false;}
    if (typeid(water_placer) != typeid(rhs.water_placer)) {return false;}
    if (width != rhs.width) {return false;}
    if (axes != rhs.axes) {return false;}
    return true;
}