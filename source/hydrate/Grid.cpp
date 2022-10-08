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
#include <data/Water.h>
#include <hydrate/AxesPlacement.h>
#include <hydrate/RadialPlacement.h>
#include <hydrate/JanPlacement.h>
#include <hydrate/CounterCulling.h>
#include <hydrate/OutlierCulling.h>
#include <hydrate/RandomCulling.h>
#include <hydrate/ClusterCulling.h>
#include <utility/Settings.h>
#include <math/Vector3.h>
#include <utility/Utility.h>

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace grid;

Grid::Grid(const Axis3D& axes, double width, double ra, double rh, setting::grid::PlacementStrategy ps, setting::grid::CullingStrategy cs) : axes(axes) {
    setup(width, ra, rh, ps, cs);
}

Grid::Grid(const vector<Atom>& atoms, double width, double ra, double rh, setting::grid::PlacementStrategy ps, setting::grid::CullingStrategy cs) {
    // find the bounding box
    Vector3<double> nmin, nmax; // new min & max
    auto[min, max] = bounding_box(atoms);

    // expand bounding box by scaling factor
    for (unsigned int i = 0; i < 3; i++) {
        double expand = 0.5*(max[i] - min[i])*setting::grid::scaling;   // amount to expand in each direction
        nmin[i] = min[i] - expand - setting::grid::width;               // ensure at least one additional bin is added
        nmax[i] = max[i] + expand + setting::grid::width;               
    }

    // setup the rest of the class members
    axes = Axis3D(nmin, nmax, setting::grid::width);
    setup(width, ra, rh, ps, cs);

    // finally add the atoms to the grid
    add(atoms);
}

Grid::Grid(const vector<Body>& bodies, double width, double ra, double rh, setting::grid::PlacementStrategy ps, setting::grid::CullingStrategy cs) {
    Vector3<double> nmin, nmax; // new min & max

    // find the total bounding box containing all bodies
    Vector3 min{0, 0, 0}, max{0, 0, 0};
    for (const Body& body : bodies) {
        auto[cmin, cmax] = bounding_box(body.get_protein_atoms());

        for (unsigned int i = 0; i < 3; i++) {
            if (cmin[i] < min[i]) {min[i] = cmin[i];}
            if (cmax[i] > max[i]) {max[i] = cmax[i];}
        }
    }

    // expand bounding box by scaling factor
    for (unsigned int i = 0; i < 3; i++) {
        double expand = 0.5*(max[i] - min[i])*setting::grid::scaling;   // amount to expand in each direction
        nmin[i] = min[i] - expand - setting::grid::width;               // ensure at least one additional bin is added
        nmax[i] = max[i] + expand + setting::grid::width;               
    }

    // setup the rest of the class members
    axes = Axis3D(nmin, nmax, setting::grid::width);
    setup(width, ra, rh, ps, cs);

    // finally add all atoms to the grid
    for (const Body& body : bodies) {
        add(body.get_protein_atoms());
    }
}

Grid::Grid(const Grid& grid) : Grid(grid.axes, grid.width, grid.ra, grid.rh, setting::grid::placement_strategy, setting::grid::culling_strategy) {
    this->grid = grid.grid;
    this->volume = grid.volume;
    this->ra = grid.ra;
    this->rh = grid.rh;
    this->a_members = grid.a_members;
    this->w_members = grid.w_members;
}

void Grid::setup(double width, double ra, double rh, setting::grid::PlacementStrategy ps, setting::grid::CullingStrategy cs) {
    // check if the grid should be cubic
    if (setting::grid::cubic) {
        double x_side = axes.x.max - axes.x.min;
        double y_side = axes.y.max - axes.y.min;
        double z_side = axes.z.max - axes.z.min;
        
        if (x_side > y_side && x_side > z_side) {
            axes.y = axes.x;
            axes.z = axes.x;
        } else if (y_side > x_side && y_side > z_side) {
            axes.x = axes.y;
            axes.z = axes.y;
        } else if (z_side > x_side && z_side > y_side) {
            axes.x = axes.z;
            axes.y = axes.z;
        }

        // update the number of bins to reflect the changed axes
        axes.rebin(width);
    }

    // check if the grid is abnormally large
    long long int total_bins = (long long) axes.x.bins*axes.y.bins*axes.z.bins;
    if (total_bins > 32e9) {
        throw except::invalid_argument("Error in Grid: Too many bins.");
    } else if (total_bins > 4e9) {
        utility::print_warning("Warning in Grid::setup: Consider lowering the number of bins.");
    }

    this->width = width;
    this->grid = GridObj(axes.x.bins, axes.y.bins, axes.z.bins);
    this->set_radius_atoms(ra);
    this->set_radius_water(rh);

    switch (ps) {
        case setting::grid::PlacementStrategy::AxesStrategy: 
            water_placer = std::make_unique<AxesPlacement>(this);
            break;
        case setting::grid::PlacementStrategy::RadialStrategy:
            water_placer = std::make_unique<RadialPlacement>(this);
            break;
        case setting::grid::PlacementStrategy::JanStrategy: 
            water_placer = std::make_unique<JanPlacement>(this);
            break;
        default: 
            throw except::unknown_argument("Error in Grid::Grid: Unkown PlacementStrategy");
    }

    switch (cs) {
        case setting::grid::CullingStrategy::CounterStrategy: 
            water_culler = std::make_unique<CounterCulling>(this);
            break;
        case setting::grid::CullingStrategy::OutlierStrategy: 
            water_culler = std::make_unique<OutlierCulling>(this);
            break;
        case setting::grid::CullingStrategy::RandomStrategy:
            water_culler = std::make_unique<RandomCulling>(this);
            break;
        default: 
            throw except::unknown_argument("Error in Grid::Grid: Unkown CullingStrategy");
    }
}

vector<Water> Grid::hydrate() {
    vector<GridMember<Water>> placed_water = find_free_locs(); // the molecules which were placed by the find_free_locs method
    water_culler->set_target_count(setting::grid::percent_water*a_members.size()); // target is 10% of atoms
    return water_culler->cull(placed_water);
}

vector<GridMember<Water>> Grid::find_free_locs() {
    // a quick check to verify there are no water molecules already present
    if (w_members.size() != 0) {
            utility::print_warning("Warning in Grid::find_free_locs: Attempting to hydrate a grid which already contains water!");
    }
    expand_volume();

    // place the water molecules with the chosen strategy
    return water_placer->place();
}

std::pair<Vector3<int>, Vector3<int>> Grid::bounding_box_index() const {
    if (__builtin_expect(a_members.size() == 0, false)) {
        throw except::invalid_operation("Error in Grid::bounding_box: Calculating a boundary box for a grid with no members!");
    }

    // initialize the bounds as large as possible
    Vector3<int> min(axes.x.bins, axes.y.bins, axes.z.bins);
    Vector3<int> max(0, 0, 0);
    for (const auto& atom : a_members) {
        for (unsigned int i = 0; i < 3; i++) {
            if (min[i] > atom.loc[i]) min[i] = atom.loc[i]; // min
            if (max[i] < atom.loc[i]) max[i] = atom.loc[i]+1; // max. +1 since this will often be used as loop limits
        }
    }
    return std::make_pair(min, max);
}

std::pair<Vector3<double>, Vector3<double>> Grid::bounding_box(const vector<Atom>& atoms) {
    // initialize the bounds as large as possible
    Vector3 min = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vector3 max = {std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min()};
    for (const auto& atom : atoms) {
        for (int i = 0; i < 3; i++) {
            min[i] = std::min(min[i], atom.coords[i]);
            max[i] = std::max(max[i], atom.coords[i]);
        }
    }
    return std::make_pair(min, max);
}

void Grid::set_radius_atoms(double radius) {
    unsigned int new_r = int(radius/width); // convert the radius to a "bin-radius"
    if (this->ra != 0 && this->ra != new_r) {
        utility::print_warning("Warning in Grid::set_radius: The radius is already set for this grid!");
    }
    this->ra = new_r;
}

void Grid::set_radius_water(double radius) {
    unsigned int new_r = int(radius/width); // convert the radius to a "bin-radius"
    if (this->rh != 0 && this->rh != new_r) {
        utility::print_warning("Warning in Grid::set_radius: The radius is already set for this grid!");
    }
    this->rh = new_r;
}

vector<Water> Grid::get_waters() const {
    vector<Water> atoms(w_members.size());
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

void Grid::expand_volume(GridMember<Water>& water) {
    if (water.expanded_volume) {return;} // check if this location has already been expanded

    water.expanded_volume = true; // mark this location as expanded
    expand_volume(water.loc, true); // do the expansion
}

void Grid::expand_volume(const Vector3<int>& loc, bool is_water) {
    char marker = is_water ? GridObj::H_AREA : GridObj::A_AREA;
    int x = loc.x(), y = loc.y(), z = loc.z();

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int r = is_water ? rh : ra; // determine which radius to use for the expansion
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) axes.x.bins); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) axes.y.bins); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) axes.z.bins); // zminus and zplus

    // loop over each bin in the box
    int added_volume = 0;

    // i, j, k *must* be ints
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(x - i, 2) + std::pow(y - j, 2) + std::pow(z - k, 2)) <= r) {
                    if (grid.index(i, j, k) != 0) {continue;} // skip if the bin is already occupied
                    added_volume++;
                    grid.index(i, j, k) = marker;
                }
            }
        }
    }

    if (!is_water) {volume += added_volume;};
}

std::vector<bool> Grid::remove_disconnected_atoms(unsigned int min) {
    expand_volume();
    ClusterCulling culler(this);
    auto to_remove = culler.remove_clusters(min);
    remove(to_remove);
    // auto to_remove_t = culler.remove_tendrils(10);
    // remove(to_remove_t);

    // // combine the two vectors
    // unsigned int index = 0;
    // for (unsigned int i = 0; i < to_remove.size(); i++) {
    //     if (!to_remove[i]) {
    //         to_remove[i] = to_remove_t[index++];
    //     }
    // }

    // // sanity check
    // if (index != to_remove_t.size()) {
    //     throw except::unexpected("Error in Grid::remove_disconnected_atoms: Index mismatch! Index is (" + std::to_string(index) + ") but the size of the vector is (" + std::to_string(to_remove_t.size()) + ")!");
    // }

    return to_remove;
}

vector<GridMember<Atom>> Grid::add(const Body* const body) {
    return add(body->get_protein_atoms());
}

void Grid::remove(const Body* const body) {
    remove(body->get_protein_atoms());
}

GridMember<Atom> Grid::add(const Atom& atom, bool expand) {
    auto loc = to_bins(atom.coords);
    unsigned int x = loc.x(), y = loc.y(), z = loc.z();

    // sanity check
    bool out_of_bounds = x >= axes.x.bins || y >= axes.y.bins || z >= axes.z.bins;
    if (__builtin_expect(out_of_bounds, false)) {
        throw except::out_of_bounds("Error in Grid::add: Atom is located outside the grid!\nBin location: " + loc.to_string() + "\n: " + axes.to_string() + "\nReal location: " + atom.coords.to_string());
    }

    if (grid.index(x, y, z) == GridObj::EMPTY) {volume++;} // can probably be removed

    GridMember gm(atom, loc);
    if (expand) {expand_volume(gm);}
    a_members.push_back(gm);
    grid.index(x, y, z) = GridObj::A_CENTER;

    return gm;
}

GridMember<Water> Grid::add(const Water& water, bool expand) {
    auto loc = to_bins(water.coords);
    unsigned int x = loc.x(), y = loc.y(), z = loc.z(); 

    // sanity check
    bool out_of_bounds = x >= axes.x.bins || y >= axes.y.bins || z >= axes.z.bins;
    if (__builtin_expect(out_of_bounds, false)) {
        throw except::out_of_bounds("Error in Grid::add: Atom is located outside the grid!\nBin location: " + loc.to_string() + "\n: " + axes.to_string() + "\nReal location: " + water.coords.to_string());
    }

    GridMember gm(water, loc);
    if (expand) {expand_volume(gm);}
    w_members.push_back(gm);
    grid.index(x, y, z) = GridObj::H_CENTER;

    return gm;
}

void Grid::remove(std::vector<bool>& to_remove) {
    // create a map of uids to remove
    std::unordered_map<unsigned int, bool> removed;

    // fill it based on the to_remove vector
    unsigned int index = 0;
    unsigned int total_removed = 0;
    for(auto& atom : a_members) {
        if (to_remove[index++]) {
            removed[atom.atom.uid] = true;
            total_removed++;
        }
    }

    index = 0;
    vector<GridMember<Atom>> removed_atoms(total_removed);
    auto predicate = [&removed, &removed_atoms, &index] (const GridMember<Atom>& gm) {
        if (removed[gm.atom.uid]) { // now we can simply look up in our removed vector to determine if an element should be removed
            removed_atoms[index++] = std::move(gm);
            return true;
        }
        return false;
    };

    // we save the sizes so we can make a sanity check after the removal    
    size_t prev_size = a_members.size();
    a_members.remove_if(predicate);
    size_t cur_size = a_members.size();

    // sanity check
    if (__builtin_expect(prev_size - cur_size != total_removed, false)) {
        throw except::invalid_operation("Error in Grid::remove: Something went wrong.");
    }

    // clean up the grid
    for (auto& atom : removed_atoms) {
        deflate_volume(atom);
        grid.index(atom.loc) = GridObj::EMPTY;
        volume--;
    }
}

void Grid::remove(const Atom& atom) {
    auto pos = std::find(a_members.begin(), a_members.end(), atom);
    if (__builtin_expect(pos == a_members.end(), false)) {
        throw except::invalid_operation("Error in Grid::remove: Attempting to remove an atom which is not part of the grid!");
    }

    GridMember<Atom>& member = *pos;
    const int x = member.loc.x(), y = member.loc.y(), z = member.loc.z();

    deflate_volume(member);
    a_members.erase(pos);
    grid.index(x, y, z) = GridObj::EMPTY;
    volume--;
}

void Grid::remove(const Water& water) {
    auto pos = std::find(w_members.begin(), w_members.end(), water);
    if (__builtin_expect(pos == w_members.end(), false)) {
        throw except::invalid_operation("Error in Grid::remove: Attempting to remove an atom which is not part of the grid!");
    }

    GridMember<Water>& member = *pos;
    const int x = member.loc.x(), y = member.loc.y(), z = member.loc.z();

    deflate_volume(member);
    w_members.erase(pos);
    grid.index(x, y, z) = GridObj::EMPTY;
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
        deflate_volume(atom);
        grid.index(atom.loc) = GridObj::EMPTY;
        volume--;
    }
}

void Grid::remove(const vector<Water>& waters) {
    // we make a vector of all possible uids
    vector<bool> removed(Atom::uid_counter);
    // and fill it with the uids that should be removed
    std::for_each(waters.begin(), waters.end(), [&removed] (const Water& water) {removed[water.uid] = true;});

    size_t index = 0; // current index in removed_waters
    vector<GridMember<Water>> removed_waters(waters.size()); // the waters which will be removed
    auto predicate = [&removed, &removed_waters, &index] (const GridMember<Water>& gm) {
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
        deflate_volume(atom);
        grid.index(atom.loc) = GridObj::EMPTY;
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

void Grid::deflate_volume(GridMember<Water>& water) {
    if (!water.expanded_volume) {return;} // check if this location has already been expanded

    water.expanded_volume = false; // mark this location as expanded
    deflate_volume(water.loc, true); // do the expansion
}

void Grid::deflate_volume(const Vector3<int>& loc, const bool is_water) {
    int x = loc.x(), y = loc.y(), z = loc.z();

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int r = is_water ? rh : ra; // determine which radius to use for the expansion
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) axes.x.bins); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) axes.y.bins); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) axes.z.bins); // zminus and zplus

    // loop over each bin in the box
    int removed_volume = -1; // -1 because we overwrite the center

    // i, j, k *must* be ints
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(loc.x() - i, 2) + std::pow(loc.y() - j, 2) + std::pow(loc.z() - k, 2)) <= r) {
                    if (grid.index(i, j, k) == GridObj::EMPTY) {continue;} // skip if bin is empty
                    removed_volume++;
                    grid.index(i, j, k) = GridObj::EMPTY;
                }
            }
        }
    }

    if (is_water) {
        grid.index(x, y, z) = GridObj::H_CENTER; // replace the center
    } else {
        grid.index(x, y, z) = GridObj::A_CENTER;
        volume -= removed_volume; // only the actual atoms contributes to the volume
    }
}

void Grid::clear_waters() {
    vector<Water> waters;
    waters.reserve(w_members.size());

    std::for_each(w_members.begin(), w_members.end(), [&waters] (const GridMember<Water>& water) {waters.push_back(water.atom);});
    remove(waters);
    assert(w_members.size() == 0);
}

Vector3<int> Grid::get_bins() const {
    return Vector3<int>(axes.x.bins, axes.y.bins, axes.z.bins);
}

Vector3<int> Grid::to_bins(const Vector3<double>& v) const {
    int binx = std::round((v.x() - axes.x.min)/width);
    int biny = std::round((v.y() - axes.y.min)/width);
    int binz = std::round((v.z() - axes.z.min)/width);
    return Vector3<int>(binx, biny, binz);
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