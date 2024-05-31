/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/detail/AtomCollection.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/GridSettings.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <constants/Constants.h>
#include <io/ExistingFile.h>

using namespace grid;
using namespace data;
using namespace data::record;

Grid::Grid(const Limit3D& axes) : axes(Axis3D(axes, settings::grid::width)) {
    setup();
}

Grid::Grid(const std::vector<Atom>& atoms) : Grid({Body(atoms)}) {}

Grid::Grid(const std::vector<Body>& bodies) {
    // find the total bounding box containing all bodies
    Vector3 min{0, 0, 0}, max{0, 0, 0};
    for (const Body& body : bodies) {
        auto[cmin, cmax] = bounding_box(body.get_atoms());

        for (unsigned int i = 0; i < 3; i++) {
            if (cmin[i] < min[i]) {min[i] = cmin[i];}
            if (cmax[i] > max[i]) {max[i] = cmax[i];}
        }
    }

    // for small systems, expand the grid by a factor 2
    auto diff = max - min;
    if (settings::grid::scaling == 0.25 && (diff.x() < 50 || diff.y() < 50 || diff.z() < 50)) {
        settings::grid::scaling = 1;
    }

    // expand bounding box by scaling factor
    Vector3<double> nmin, nmax; // new min & max
    for (unsigned int i = 0; i < 3; i++) {
        double expand = 0.5*diff[i]*settings::grid::scaling; // amount to expand in each direction
        nmin[i] = std::floor(min[i] - expand - settings::grid::width); // flooring to make our grid 'nice' (i.e. bin edges are at integer values)
        nmax[i] = std::ceil( max[i] + expand + settings::grid::width); //  ceiling to make our grid 'nice' (i.e. bin edges are at integer values)
    }

    // setup the rest of the class members
    axes = Axis3D(nmin, nmax, settings::grid::width);
    setup();

    // finally add all atoms to the grid
    for (const Body& body : bodies) {
        add(body.get_atoms());
    }
}

Grid::Grid(const Grid& grid) {
    *this = grid;
}

Grid::Grid(Grid&& grid) noexcept {
    *this = std::move(grid);
}

Grid::~Grid() = default;

void Grid::setup() {
    // check if the grid should be cubic
    if (settings::grid::cubic) {
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
    }

    // enforce minimum number of bins if set
    if (settings::grid::min_bins != 0) {
        double min_length = settings::grid::min_bins*settings::grid::width/2;
        if (axes.x.bins < settings::grid::min_bins) {
            axes.x.bins = settings::grid::min_bins;
            axes.x.min = -min_length;
            axes.x.max =  min_length;
        }

        if (axes.y.bins < settings::grid::min_bins) {
            axes.y.bins = settings::grid::min_bins;
            axes.y.min = -min_length;
            axes.y.max =  min_length;
        }

        if (axes.z.bins < settings::grid::min_bins) {
            axes.z.bins = settings::grid::min_bins;
            axes.z.min = -min_length;
            axes.z.max =  min_length;
        }
    }

    // check if the grid is abnormally large
    long long int total_bins = (long long) axes.x.bins*axes.y.bins*axes.z.bins;
    if (total_bins > 32e9) {
        throw except::size_error("Grid::setup: Attempting to allocate a grid of size > 16GB. Try reducing the number of bins.");
    } else if (total_bins > 4e9) {
        console::print_warning("Warning in Grid::setup: Attempting to allocate a grid of size > 2GB. Consider lowering the number of bins.");
    }

    this->grid = detail::GridObj(axes.x.bins, axes.y.bins, axes.z.bins);
}

double Grid::get_atomic_radius(constants::atom_t atom) const {
    return constants::radius::get_vdw_radius(atom);
}

double Grid::get_hydration_radius() const {
    return constants::radius::get_vdw_radius(constants::atom_t::O);
}

std::pair<Vector3<int>, Vector3<int>> Grid::bounding_box_index() const {
    if (a_members.size() == 0) [[unlikely]] {
        throw except::invalid_operation("Grid::bounding_box: Calculating a boundary box for a grid with no members!");
    }

    // initialize the bounds as large as possible
    Vector3<int> min(axes.x.bins, axes.y.bins, axes.z.bins);
    Vector3<int> max(0, 0, 0);
    for (const auto& atom : a_members) {
        for (unsigned int i = 0; i < 3; i++) {
            if (min[i] > atom.get_bin_loc()[i]) min[i] = atom.get_bin_loc()[i]; // min
            if (max[i] < atom.get_bin_loc()[i]) max[i] = atom.get_bin_loc()[i]+1; // max. +1 since this will often be used as loop limits
        }
    }
    return std::make_pair(min, max);
}

std::pair<Vector3<double>, Vector3<double>> Grid::bounding_box(const std::vector<Atom>& atoms) {
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

std::vector<Water> Grid::get_waters() const {
    std::vector<Water> atoms(w_members.size());
    int i = 0; // counter
    for (const auto& water : w_members) {
        atoms[i] = water.get_atom();
        i++;
    }
    atoms.resize(i);
    return atoms;
}

std::vector<Atom> Grid::get_atoms() const {
    std::vector<Atom> atoms(a_members.size());
    int i = 0; // counter
    for (const auto& atom : a_members) {
        atoms[i] = atom.get_atom();
        i++;
    }
    atoms.resize(i);
    return atoms;
}

void Grid::add_volume(double value) {
    volume += value;
}

void Grid::force_expand_volume() {
    for (auto& atom : a_members) {
        atom.set_expanded(false);
        expand_volume(atom);
    }

    for (auto& water : w_members) {
        water.set_expanded(false);
        expand_volume(water);
    }
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
    if (atom.is_expanded()) {return;} // check if this location has already been expanded
    atom.set_expanded(true); // mark this location as expanded

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = atom.get_bin_loc().x(), y = atom.get_bin_loc().y(), z = atom.get_bin_loc().z(); 
    double rvdw = get_atomic_radius(atom.get_atom_type())/settings::grid::width;
    double rvol = settings::grid::rvol/settings::grid::width;
    double rvdw2 = std::pow(rvdw, 2), rvol2 = std::pow(rvol, 2);

    int xm = std::max<int>(x - std::ceil(rvol), 0), xp = std::min<int>(x + std::ceil(rvol) + 1, axes.x.bins); // xminus and xplus
    int ym = std::max<int>(y - std::ceil(rvol), 0), yp = std::min<int>(y + std::ceil(rvol) + 1, axes.y.bins); // yminus and yplus
    int zm = std::max<int>(z - std::ceil(rvol), 0), zp = std::min<int>(z + std::ceil(rvol) + 1, axes.z.bins); // zminus and zplus

    // loop over each bin in the box
    int added_volume = 0;

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(x - i, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(y - j, 2);
            for (int k = zm; k < zp; ++k) {
                // fill a sphere of radius [0, vdw] around the atom
                double dist = x2y2 + std::pow(z - k, 2);
                auto& bin = grid.index(i, j, k);
                if (dist <= rvdw2) {
                    if (!grid.is_empty_or_volume(bin)) {continue;}
                    added_volume += !grid.is_volume(bin); // only add to the volume if the bin is not already part of the volume
                    bin = detail::A_AREA;
                }

                // fill an outer shell of radius [vdw, rvol] to make sure the volume is space-filling
                else if (dist <= rvol2) {
                    if (!grid.is_empty(i, j, k)) {continue;}
                    bin = detail::VOLUME;
                    added_volume++;
                }
            }
        }
    }
    volume += added_volume;
}

void Grid::expand_volume(GridMember<Water>& water) {
    if (water.is_expanded()) {return;} // check if this location has already been expanded
    water.set_expanded(true); // mark this location as expanded

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = water.get_bin_loc().x(), y = water.get_bin_loc().y(), z = water.get_bin_loc().z();
    double rvdw = get_hydration_radius()/settings::grid::width;
    double rvdw2 = std::pow(rvdw, 2);
    int xm = std::max<int>(x - std::ceil(rvdw), 0), xp = std::min<int>(x + std::ceil(rvdw) + 1, axes.x.bins); // xminus and xplus
    int ym = std::max<int>(y - std::ceil(rvdw), 0), yp = std::min<int>(y + std::ceil(rvdw) + 1, axes.y.bins); // yminus and yplus
    int zm = std::max<int>(z - std::ceil(rvdw), 0), zp = std::min<int>(z + std::ceil(rvdw) + 1, axes.z.bins); // zminus and zplus

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(x - i, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(y - j, 2);
            for (int k = zm; k < zp; ++k) {
                // determine if the bin is within a sphere centered on the atom
                double dist = x2y2 + std::pow(z - k, 2);
                if (dist <= rvdw2) {
                    if (!grid.is_empty(i, j, k)) {continue;}
                    grid.index(i, j, k) = detail::W_AREA;
                }
            }
        }
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
    if (!atom.is_expanded()) {return;} // check if this location has already been deflated
    atom.set_expanded(false); // mark the atom as deflated

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = atom.get_bin_loc().x(), y = atom.get_bin_loc().y(), z = atom.get_bin_loc().z();
    double rvol = std::max<double>(get_atomic_radius(atom.get_atom_type()), settings::grid::rvol)/settings::grid::width;
    double rvol2 = std::pow(rvol, 2);
    int xm = std::max<int>(x - std::ceil(rvol), 0), xp = std::min<int>(x + std::ceil(rvol) + 1, axes.x.bins); // xminus and xplus
    int ym = std::max<int>(y - std::ceil(rvol), 0), yp = std::min<int>(y + std::ceil(rvol) + 1, axes.y.bins); // yminus and yplus
    int zm = std::max<int>(z - std::ceil(rvol), 0), zp = std::min<int>(z + std::ceil(rvol) + 1, axes.z.bins); // zminus and zplus

    // i, j, k *must* be ints due to avoid unsigned underflow
    int removed_volume = 0;
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(x - i, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(y - j, 2);
            for (int k = zm; k < zp; ++k) {
                double dist = x2y2 + std::pow(z - k, 2);

                // determine if the bin is within a sphere centered on the atom
                auto& bin = grid.index(i, j, k);
                if (dist <= rvol2) {
                    if (!grid.is_atom_area_or_volume(bin)) {continue;}
                    bin = detail::EMPTY;
                    removed_volume++;
                }
            }
        }
    }
    volume -= removed_volume; // only the actual atoms contributes to the volume
}

void Grid::deflate_volume(GridMember<Water>& water) {
    if (!water.is_expanded()) {return;} // check if this location has already been deflated
    water.set_expanded(false); // mark the water as deflated

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = water.get_bin_loc().x(), y = water.get_bin_loc().y(), z = water.get_bin_loc().z();
    double rvdw = get_hydration_radius()/settings::grid::width;
    double rvdw2 = std::pow(rvdw, 2);
    int xm = std::max<int>(x - std::ceil(rvdw), 0), xp = std::min<int>(x + std::ceil(rvdw) + 1, axes.x.bins); // xminus and xplus
    int ym = std::max<int>(y - std::ceil(rvdw), 0), yp = std::min<int>(y + std::ceil(rvdw) + 1, axes.y.bins); // yminus and yplus
    int zm = std::max<int>(z - std::ceil(rvdw), 0), zp = std::min<int>(z + std::ceil(rvdw) + 1, axes.z.bins); // zminus and zplus

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(x - i, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(y - j, 2);
            for (int k = zm; k < zp; ++k) {
                double dist = x2y2 + std::pow(z - k, 2);

                // determine if the bin is within a sphere centered on the atom
                if (dist <= rvdw2) {
                    if (!grid.is_water_area(i, j, k)) {continue;}
                    grid.index(i, j, k) = detail::EMPTY;
                }
            }
        }
    }
}

#include <hydrate/culling/ClusterCulling.h>
std::vector<bool> Grid::remove_disconnected_atoms(unsigned int min) {
    // throw except::not_implemented("Grid::remove_disconnected_atoms: Not implemented!");
    expand_volume();
    ClusterCulling culler(this);
    auto to_remove = culler.remove_clusters(min);
    remove(to_remove);

    auto to_remove_t = culler.remove_tendrils(10);
    remove(to_remove_t);

    // combine the two vectors
    unsigned int index = 0;
    for (unsigned int i = 0; i < to_remove.size(); i++) {
        if (!to_remove[i]) {
            to_remove[i] = to_remove_t[index++];
        }
    }

    // sanity check
    if (index != to_remove_t.size()) {
        throw except::unexpected("Grid::remove_disconnected_atoms: Index mismatch! Index is (" + std::to_string(index) + ") but the size of the vector is (" + std::to_string(to_remove_t.size()) + ")!");
    }

    return to_remove;
}

template <typename T, typename>
std::vector<GridMember<T>> Grid::add(const std::vector<T>& atoms) {
    std::vector<GridMember<T>> v(atoms.size());
    unsigned int index = 0;
    for (const auto& a : atoms) {
        v[index++] = add(a);
    }
    return v;
}
template std::vector<GridMember<Atom>> Grid::add<Atom>(const std::vector<Atom>&);
template std::vector<GridMember<Water>> Grid::add<Water>(const std::vector<Water>&);

std::vector<GridMember<Atom>> Grid::add(const Body* body) {
    return add(body->get_atoms());
}

void Grid::remove(const Body* body) {
    remove(body->get_atoms());
}

const GridMember<Atom>& Grid::add(const Atom& atom, bool expand) {
    auto loc = to_bins(atom.coords);
    unsigned int i = loc.x(), j = loc.y(), k = loc.z();

    // sanity check
    #if DEBUG
        bool out_of_bounds = i >= axes.x.bins || j >= axes.y.bins || k >= axes.z.bins;
        if (out_of_bounds) [[unlikely]] {
            throw except::out_of_bounds("Grid::add: Atom is located outside the grid!\nBin location: " + loc.to_string() + "\n: " + axes.to_string() + "\nReal location: " + atom.coords.to_string());
        }
    #endif

    auto& bin = grid.index(i, j, k);
    volume += grid.is_empty(bin);
    bin = detail::A_CENTER;

    GridMember gm(atom, loc);
    if (expand) {expand_volume(gm);}
    a_members.push_back(std::move(gm));

    return a_members.back();
}

const GridMember<Water>& Grid::add(const Water& water, bool expand) {
    auto loc = to_bins(water.coords);
    int i = loc.x(), j = loc.y(), k = loc.z(); 

    // sanity check
    #if DEBUG
        if (i >= (int) axes.x.bins || j >= (int) axes.y.bins || k >= (int) axes.z.bins || i < 0 || j < 0 || k < 0) [[unlikely]] {
            throw except::out_of_bounds("Grid::add: Atom is located outside the grid!\nBin location: " + loc.to_string() + "\n: " + axes.to_string() + "\nReal location: " + water.coords.to_string());
        } 

        if (!(grid.index(i, j, k) == detail::EMPTY || grid.index(i, j, k) == detail::W_CENTER || grid.index(i, j, k) == detail::VOLUME)) {
            throw except::invalid_operation("Grid::add: Attempting to add a water molecule to a non-empty location (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
        }
    #endif

    GridMember gm(water, loc);
    grid.index(i, j, k) = detail::W_CENTER;
    if (expand) {expand_volume(gm);}
    w_members.push_back(std::move(gm));

    return w_members.back();
}

void Grid::remove(std::vector<bool>& to_remove) {
    // create a map of uids to remove
    std::unordered_map<unsigned int, bool> removed;

    // fill it based on the to_remove vector
    unsigned int index = 0;
    unsigned int total_removed = 0;
    for(auto& atom : a_members) {
        if (to_remove[index++]) {
            removed[atom.get_atom().uid] = true;
            total_removed++;
        } else {
            removed[atom.get_atom().uid] = false;
        }
    }

    index = 0;
    std::vector<GridMember<Atom>> removed_atoms(total_removed);
    auto predicate = [&removed, &removed_atoms, &index] (GridMember<Atom>& gm) {
        if (removed[gm.get_atom().uid]) { // now we can simply look up in our removed vector to determine if an element should be removed
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
    if (prev_size - cur_size != total_removed) [[unlikely]] {
        throw except::invalid_operation("Grid::remove: Something went wrong.");
    }

    // clean up the grid
    for (auto& atom : removed_atoms) {
        deflate_volume(atom);
        grid.index(atom.get_bin_loc()) = detail::EMPTY;
        volume--;
    }
}

void Grid::remove(const Atom& atom) {
    auto pos = std::find(a_members.begin(), a_members.end(), atom);
    if (pos == a_members.end()) [[unlikely]] {
        throw except::invalid_operation("Grid::remove: Attempting to remove an atom which is not part of the grid!");
    }

    GridMember<Atom>& member = *pos;

    a_members.erase(pos);
    deflate_volume(member);
    grid.index(member.get_bin_loc()) = detail::EMPTY;
    volume--;
}

void Grid::remove(const Water& water) {
    auto pos = std::find(w_members.begin(), w_members.end(), water);
    if (pos == w_members.end()) [[unlikely]] {
        throw except::invalid_operation("Grid::remove: Attempting to remove an atom which is not part of the grid!");
    }

    GridMember<Water>& member = *pos;

    w_members.erase(pos);
    deflate_volume(member);
    grid.index(member.get_bin_loc()) = detail::EMPTY;
}

void Grid::remove(const std::vector<Atom>& atoms) {
    // we make a vector of all possible uids
    std::unordered_map<int, bool> removed;
    // and fill it with the uids that should be removed
    std::for_each(atoms.begin(), atoms.end(), [&removed] (const Atom& atom) {removed[atom.uid] = true;});

    unsigned int index = 0; // current index in removed_atoms
    std::vector<GridMember<Atom>> removed_atoms(atoms.size()); // the atoms which will be removed
    auto predicate = [&removed, &removed_atoms, &index] (const GridMember<Atom>& gm) {
        if (removed[gm.get_atom().uid]) { // now we can simply look up in our removed vector to determine if an element should be removed
            removed_atoms[index++] = gm;
            return true;
        }
        return false;
    };

    // we save the sizes so we can make a sanity check after the removal    
    unsigned int prev_size = a_members.size();
    a_members.remove_if(predicate);
    unsigned int cur_size = a_members.size();

    // sanity check
    if (prev_size - cur_size != atoms.size()) [[unlikely]] {
        throw except::unexpected("Grid::remove: Expected to remove " + std::to_string(atoms.size()) + " elements, but only " + std::to_string(prev_size - cur_size) + " were actually removed.");
    }

    // clean up the grid
    for (auto& atom : removed_atoms) {
        deflate_volume(atom);
        auto& bin = grid.index(atom.get_bin_loc());
        volume -= grid.is_atom_center(bin); // bin may in rare cases contain two atoms, so we need to check
        bin = detail::EMPTY;
    }
}

void Grid::remove(const std::vector<Water>& waters) {
    // we make a map and fill it with the uids that should be removed
    std::unordered_map<unsigned int, bool> to_remove_id;
    std::for_each(waters.begin(), waters.end(), [&to_remove_id] (const Water& water) {to_remove_id[water.uid] = true;});

    unsigned int index = 0; // current index in removed_waters
    std::vector<GridMember<Water>> to_remove(waters.size()); // the waters which will be removed
    auto predicate = [&to_remove_id, &to_remove, &index] (const GridMember<Water>& gm) {
        if (to_remove_id.contains(gm.get_atom().uid)) { // now we can simply look up in our removed vector to determine if an element should be removed
            to_remove[index++] = gm;
            return true;
        }
        return false;
    };

    // we save the sizes so we can make a sanity check after the removal    
    unsigned int prev_size = w_members.size();
    w_members.remove_if(predicate);
    unsigned int cur_size = w_members.size();

    // sanity check
    if (prev_size - cur_size != waters.size()) [[unlikely]] {
        throw except::unexpected("Grid::remove: Expected to remove " + std::to_string(waters.size()) + " elements, but only " + std::to_string(prev_size - cur_size) + " were actually removed.");
    }

    // clean up the grid
    for (auto& water : to_remove) {
        deflate_volume(water);
        grid.index(water.get_bin_loc()) = detail::EMPTY;
    }
}

void Grid::clear_waters() {
    std::vector<Water> waters;
    waters.reserve(w_members.size());
    std::for_each(w_members.begin(), w_members.end(), [&waters] (const GridMember<Water>& water) {waters.push_back(water.get_atom());});
    remove(waters);
    if (w_members.size() != 0) [[unlikely]] {throw except::unexpected("Grid::clear_waters: Something went wrong.");}
}

Vector3<int> Grid::get_bins() const {
    return Vector3<int>(axes.x.bins, axes.y.bins, axes.z.bins);
}

Vector3<int> Grid::to_bins(const Vector3<double>& v) const {
    int binx = std::round((v.x() - axes.x.min)/settings::grid::width);
    int biny = std::round((v.y() - axes.y.min)/settings::grid::width);
    int binz = std::round((v.z() - axes.z.min)/settings::grid::width);
    return Vector3<int>(binx, biny, binz);
}

Vector3<int> Grid::to_bins_bounded(const Vector3<double>& v) const {
    auto bins = to_bins(v);
    bins.x() = std::clamp<int>(bins.x(), 0, axes.x.bins-1);
    bins.y() = std::clamp<int>(bins.y(), 0, axes.y.bins-1);
    bins.z() = std::clamp<int>(bins.z(), 0, axes.z.bins-1);
    return bins;
}

double Grid::get_volume() {
    expand_volume();
    return volume*std::pow(settings::grid::width, 3);
}

Grid& Grid::operator=(const Grid& rhs) {
    grid = rhs.grid;
    a_members = rhs.a_members;
    w_members = rhs.w_members;
    volume = rhs.volume;
    axes = rhs.axes;
    // culler & placer cannot be modified after program is run, so they'll automatically be equal always
    return *this;
}

Grid& Grid::operator=(Grid&& rhs) noexcept {
    grid = std::move(rhs.grid);
    a_members = std::move(rhs.a_members);
    w_members = std::move(rhs.w_members);
    volume = rhs.volume;
    axes = std::move(rhs.axes);
    return *this;
}

bool Grid::operator==(const Grid& rhs) const {
    // we do everything but check the contents of the grid. 
    if (volume != rhs.volume) {return false;}
    if (a_members.size() != rhs.a_members.size()) {return false;}
    if (w_members.size() != rhs.w_members.size()) {return false;}
    if (axes != rhs.axes) {return false;}
    return true;
}

void Grid::save(const io::File& path) const {
    std::vector<Atom> atoms;
    std::vector<Water> waters;
    unsigned int c = 0;
    for (unsigned int i = 0; i < grid.size_x(); i++) {
        for (unsigned int j = 0; j < grid.size_y(); j++) {
            for (unsigned int k = 0; k < grid.size_z(); k++) {
                switch (grid.index(i, j, k)) {
                    case detail::A_CENTER:
                        atoms.push_back(Atom(c++, "C", "", "LYS", 'A', 1, "", to_xyz(i, j, k), 1, 0, constants::atom_t::C, ""));
                        break;
                    case detail::A_AREA:
                        atoms.push_back(Atom(c++, "C", "", "LYS", 'B', 2, "", to_xyz(i, j, k), 1, 0, constants::atom_t::C, ""));
                        break;
                    case detail::VOLUME:
                        atoms.push_back(Atom(c++, "C", "", "LYS", 'C', 3, "", to_xyz(i, j, k), 1, 0, constants::atom_t::C, ""));
                        break;
                    case detail::W_CENTER:
                        waters.push_back(Atom(c++, "C", "", "LYS", 'D', 4, "", to_xyz(i, j, k), 1, 0, constants::atom_t::C, ""));
                        break;
                    case detail::W_AREA:
                        waters.push_back(Atom(c++, "C", "", "LYS", 'E', 5, "", to_xyz(i, j, k), 1, 0, constants::atom_t::C, ""));
                        break;
                    default:
                        break;
                }
            }
        }
    }

    // visualize corners
    atoms.push_back(Atom(c++, "C", "", "LYS", 'F', 6, "", to_xyz(0, 0, 0), 1, 0, constants::atom_t::C, ""));
    atoms.push_back(Atom(c++, "C", "", "LYS", 'F', 6, "", to_xyz(axes.x.bins, 0, 0), 1, 0, constants::atom_t::C, ""));
    atoms.push_back(Atom(c++, "C", "", "LYS", 'F', 6, "", to_xyz(axes.x.bins, axes.y.bins, 0), 1, 0, constants::atom_t::C, ""));
    atoms.push_back(Atom(c++, "C", "", "LYS", 'F', 6, "", to_xyz(0, axes.y.bins, 0), 1, 0, constants::atom_t::C, ""));
    atoms.push_back(Atom(c++, "C", "", "LYS", 'F', 6, "", to_xyz(0, axes.y.bins, axes.z.bins), 1, 0, constants::atom_t::C, ""));
    atoms.push_back(Atom(c++, "C", "", "LYS", 'F', 6, "", to_xyz(0, 0, axes.z.bins), 1, 0, constants::atom_t::C, ""));
    atoms.push_back(Atom(c++, "C", "", "LYS", 'F', 6, "", to_xyz(axes.x.bins, 0, axes.z.bins), 1, 0, constants::atom_t::C, ""));
    atoms.push_back(Atom(c++, "C", "", "LYS", 'F', 6, "", to_xyz(axes.x.bins, axes.y.bins, axes.z.bins), 1, 0, constants::atom_t::C, ""));

    data::Molecule p(atoms, waters);
    p.save(path);
}

grid::detail::GridExcludedVolume Grid::generate_excluded_volume(bool determine_surface) {
    expand_volume();
    std::vector<Vector3<double>> exv_atoms;
    exv_atoms.reserve(volume);
    auto[vmin, vmax] = bounding_box_index();

    if (std::abs(std::remainder(settings::grid::exv_radius*2, settings::grid::width)) > 1e-6) {
        console::print_warning("Warning in Grid::generate_excluded_volume: The excluded volume radius is not a multiple of the grid width. Rounding to nearest integer.");
    }
    int stride = std::round(2*settings::grid::exv_radius/settings::grid::width);
    int buffer = 2./settings::grid::width; // 2Å buffer in each direction should be enough to capture all filled voxels
    for (int i = std::max<int>(vmin.x()-buffer, 0); i < std::min<int>(vmax.x()+buffer, axes.x.bins); i+=stride) {
        for (int j = std::max<int>(vmin.y()-buffer, 0); j < std::min<int>(vmax.y()+buffer, axes.y.bins); j+=stride) {
            for (int k = std::max<int>(vmin.z()-buffer, 0); k < std::min<int>(vmax.z()+buffer, axes.z.bins); k+=stride) {
                switch (grid.index(i, j, k)) {
                    case detail::VOLUME:
                    case detail::A_AREA:
                    case detail::A_CENTER: 
                        exv_atoms.push_back(to_xyz(i, j, k));
                    default:
                        break;
                }
            }
        }
    }

    if (settings::grid::save_exv) {
        std::vector<Atom> atoms(exv_atoms.size());
        for (unsigned int i = 0; i < exv_atoms.size(); i++) {
            atoms[i] = Atom(i, "C", "", "LYS", 'A', 1, "", exv_atoms[i], 1, 0, constants::atom_t::C, "");
        }
        data::detail::AtomCollection(atoms, {}).write(settings::general::output + "exv.pdb");
    }

    if (!determine_surface) {
        return {exv_atoms, {}};
    }

    // determine the surface atoms
    std::vector<Vector3<double>> surface_atoms;
    
}

const grid::detail::State& Grid::index(unsigned int i, unsigned int j, unsigned int k) const {
    return grid.index(i, j, k);
}

std::vector<Atom> Grid::get_surface_atoms() const {
    throw except::not_implemented("Grid::get_surface_atoms: Not implemented.");
}

Vector3<double> Grid::to_xyz(int i, int j, int k) const {
    double x = axes.x.min + i*settings::grid::width;
    double y = axes.y.min + j*settings::grid::width;
    double z = axes.z.min + k*settings::grid::width;
    return {x, y, z};
}

Vector3<int> Grid::get_center() const {
    return {int(axes.x.bins/2), int(axes.y.bins/2), int(axes.z.bins/2)};
}

double Grid::get_width() const {return settings::grid::width;}