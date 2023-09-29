#include <hydrate/Grid.h>
#include <hydrate/culling/CullingFactory.h>
#include <hydrate/placement/PlacementFactory.h>
#include <hydrate/GridMember.h>
#include <data/Atom.h>
#include <data/Body.h>
#include <data/Water.h>
#include <settings/GridSettings.h>
#include <utility/Console.h>
#include <utility/Constants.h>
#include <data/Protein.h>
#include <io/ExistingFile.h>

using namespace grid;

Grid::Grid(const Axis3D& axes) : axes(axes) {
    setup();
}

Grid::Grid(const std::vector<Atom>& atoms) : Grid({Body(atoms)}) {}

Grid::Grid(const std::vector<Body>& bodies) {
    Vector3<double> nmin, nmax; // new min & max

    // find the total bounding box containing all bodies
    Vector3 min{0, 0, 0}, max{0, 0, 0};
    for (const Body& body : bodies) {
        auto[cmin, cmax] = bounding_box(body.get_atoms());

        for (unsigned int i = 0; i < 3; i++) {
            if (cmin[i] < min[i]) {min[i] = std::floor(cmin[i]);}
            if (cmax[i] > max[i]) {max[i] = std::ceil(cmax[i]);}
        }
    }

    // expand bounding box by scaling factor
    for (unsigned int i = 0; i < 3; i++) {
        double expand = 0.5*(max[i] - min[i])*settings::grid::scaling;   // amount to expand in each direction
        nmin[i] = min[i] - expand - settings::grid::width;               // ensure at least one additional bin is added
        nmax[i] = max[i] + expand + settings::grid::width;               
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
    // make sure we have the correct number of bins
    axes.rebin(settings::grid::width);

    // check if the grid is abnormally large
    long long int total_bins = (long long) axes.x.bins*axes.y.bins*axes.z.bins;
    if (total_bins > 32e9) {
        throw except::invalid_argument("Grid: Too many bins.");
    } else if (total_bins > 4e9) {
        console::print_warning("Warning in Grid::setup: Consider lowering the number of bins.");
    }

    this->width = settings::grid::width;
    this->grid = GridObj(axes.x.bins, axes.y.bins, axes.z.bins);
    this->set_radius_atoms(settings::grid::ra);
    this->set_radius_water(settings::grid::rh);

    water_placer = grid::factory::construct_placement_strategy(this);
    water_culler = grid::factory::construct_culling_strategy(this);
}

std::vector<Water> Grid::hydrate() {
    std::vector<GridMember<Water>> placed_water = find_free_locs(); // the molecules which were placed by the find_free_locs method
    water_culler->set_target_count(settings::grid::percent_water*a_members.size()); // target is 10% of atoms
    return water_culler->cull(placed_water);
}

std::vector<GridMember<Water>> Grid::find_free_locs() {
    // a quick check to verify there are no water molecules already present
    if (w_members.size() != 0) {
            console::print_warning("Warning in Grid::find_free_locs: Attempting to hydrate a grid which already contains water!");
    }
    expand_volume();

    // place the water molecules with the chosen strategy
    return water_placer->place();
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
            if (min[i] > atom.loc[i]) min[i] = atom.loc[i]; // min
            if (max[i] < atom.loc[i]) max[i] = atom.loc[i]+1; // max. +1 since this will often be used as loop limits
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

void Grid::set_radius_atoms(double radius) {
    unsigned int new_r = int(radius/width); // convert the radius to a "bin-radius"
    if (this->ra != 0 && this->ra != new_r) {
        console::print_warning("Warning in Grid::set_radius: The radius is already set for this grid!");
    }
    this->ra = new_r;
}

void Grid::set_radius_water(double radius) {
    unsigned int new_r = int(radius/width); // convert the radius to a "bin-radius"
    if (this->rh != 0 && this->rh != new_r) {
        console::print_warning("Warning in Grid::set_radius: The radius is already set for this grid!");
    }
    this->rh = new_r;
}

std::vector<Water> Grid::get_waters() const {
    std::vector<Water> atoms(w_members.size());
    int i = 0; // counter
    for (const auto& water : w_members) {
        atoms[i] = water.atom;
        i++;
    }
    atoms.resize(i);
    return atoms;
}

std::vector<Atom> Grid::get_atoms() const {
    std::vector<Atom> atoms(a_members.size());
    int i = 0; // counter
    for (const auto& atom : a_members) {
        atoms[i] = atom.atom;
        i++;
    }
    atoms.resize(i);
    return atoms;
}

void Grid::force_expand_volume() {
    for (auto& atom : a_members) {
        atom.expanded_volume = false;
        expand_volume(atom);
    }

    for (auto& water : w_members) {
        water.expanded_volume = false;
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
    if (atom.expanded_volume) {return;} // check if this location has already been expanded
    atom.expanded_volume = true; // mark this location as expanded

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = atom.loc.x(), y = atom.loc.y(), z = atom.loc.z();
    int r = ra; // to avoid underflow
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) axes.x.bins); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) axes.y.bins); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) axes.z.bins); // zminus and zplus

    // loop over each bin in the box
    int added_volume = 0;

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(x - i, 2) + std::pow(y - j, 2) + std::pow(z - k, 2)) <= r) {
                    if (grid.index(i, j, k) != GridObj::EMPTY) {continue;}
                    grid.index(i, j, k) = GridObj::A_AREA;
                    added_volume++;
                }
            }
        }
    }
    volume += added_volume;
}

void Grid::expand_volume(GridMember<Water>& water) {
    if (water.expanded_volume) {return;} // check if this location has already been expanded
    water.expanded_volume = true; // mark this location as expanded

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = water.loc.x(), y = water.loc.y(), z = water.loc.z();
    int r = rh; // to avoid unsigned underflow
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) axes.x.bins); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) axes.y.bins); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) axes.z.bins); // zminus and zplus

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(x - i, 2) + std::pow(y - j, 2) + std::pow(z - k, 2)) <= r) {
                    if (grid.index(i, j, k) != GridObj::EMPTY) {continue;}
                    grid.index(i, j, k) = GridObj::H_AREA;
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

template <typename T, typename = std::enable_if_t<std::is_base_of<Atom, T>::value>>
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
    unsigned int x = loc.x(), y = loc.y(), z = loc.z();

    // sanity check
    bool out_of_bounds = x >= axes.x.bins || y >= axes.y.bins || z >= axes.z.bins;
    if (out_of_bounds) [[unlikely]] {
        throw except::out_of_bounds("Grid::add: Atom is located outside the grid!\nBin location: " + loc.to_string() + "\n: " + axes.to_string() + "\nReal location: " + atom.coords.to_string());
    }

    if (grid.index(x, y, z) != GridObj::A_AREA) {volume++;}

    GridMember gm(atom, loc);
    grid.index(x, y, z) = GridObj::A_CENTER;
    if (expand) {expand_volume(gm);}
    a_members.push_back(std::move(gm));

    return a_members.back();
}

const GridMember<Water>& Grid::add(const Water& water, bool expand) {
    auto loc = to_bins(water.coords);
    unsigned int x = loc.x(), y = loc.y(), z = loc.z(); 

    // sanity check
    bool out_of_bounds = x >= axes.x.bins || y >= axes.y.bins || z >= axes.z.bins;
    if (out_of_bounds) [[unlikely]] {
        throw except::out_of_bounds("Grid::add: Atom is located outside the grid!\nBin location: " + loc.to_string() + "\n: " + axes.to_string() + "\nReal location: " + water.coords.to_string());
    }

    GridMember gm(water, loc);
    grid.index(x, y, z) = GridObj::H_CENTER;
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
            removed[atom.atom.uid] = true;
            total_removed++;
        } else {
            removed[atom.atom.uid] = false;
        }
    }

    index = 0;
    std::vector<GridMember<Atom>> removed_atoms(total_removed);
    auto predicate = [&removed, &removed_atoms, &index] (GridMember<Atom>& gm) {
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
    if (prev_size - cur_size != total_removed) [[unlikely]] {
        throw except::invalid_operation("Grid::remove: Something went wrong.");
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
    if (pos == a_members.end()) [[unlikely]] {
        throw except::invalid_operation("Grid::remove: Attempting to remove an atom which is not part of the grid!");
    }

    GridMember<Atom>& member = *pos;

    a_members.erase(pos);
    deflate_volume(member);
    grid.index(member.loc) = GridObj::EMPTY;
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
    grid.index(member.loc) = GridObj::EMPTY;
}

void Grid::remove(const std::vector<Atom>& atoms) {
    // we make a vector of all possible uids
    std::unordered_map<int, bool> removed;
    // and fill it with the uids that should be removed
    std::for_each(atoms.begin(), atoms.end(), [&removed] (const Atom& atom) {removed[atom.uid] = true;});

    unsigned int index = 0; // current index in removed_atoms
    std::vector<GridMember<Atom>> removed_atoms(atoms.size()); // the atoms which will be removed
    auto predicate = [&removed, &removed_atoms, &index] (const GridMember<Atom>& gm) {
        if (removed[gm.atom.uid]) { // now we can simply look up in our removed vector to determine if an element should be removed
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
        grid.index(atom.loc) = GridObj::EMPTY;
        volume--;
    }
}

void Grid::remove(const std::vector<Water>& waters) {
    // we make a map and fill it with the uids that should be removed
    std::unordered_map<unsigned int, bool> to_remove_id;
    std::for_each(waters.begin(), waters.end(), [&to_remove_id] (const Water& water) {to_remove_id[water.uid] = true;});

    unsigned int index = 0; // current index in removed_waters
    std::vector<GridMember<Water>> to_remove(waters.size()); // the waters which will be removed
    auto predicate = [&to_remove_id, &to_remove, &index] (const GridMember<Water>& gm) {
        if (to_remove_id.contains(gm.atom.uid)) { // now we can simply look up in our removed vector to determine if an element should be removed
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
        grid.index(water.loc) = GridObj::EMPTY;
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
    if (!atom.expanded_volume) {return;} // check if this location has already been deflated
    atom.expanded_volume = false; // mark the atom as deflated

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = atom.loc.x(), y = atom.loc.y(), z = atom.loc.z();
    int r = ra; // to avoid unsigned underflow
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) axes.x.bins); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) axes.y.bins); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) axes.z.bins); // zminus and zplus

    // i, j, k *must* be ints due to avoid unsigned underflow
    int removed_volume = 0;
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(x - i, 2) + std::pow(y - j, 2) + std::pow(z - k, 2)) <= ra) {
                    if (grid.index(i, j, k) != GridObj::A_AREA) {continue;}
                    grid.index(i, j, k) = GridObj::EMPTY;
                    removed_volume++;
                }
            }
        }
    }
    volume -= removed_volume; // only the actual atoms contributes to the volume
}

void Grid::deflate_volume(GridMember<Water>& water) {
    if (!water.expanded_volume) {return;} // check if this location has already been deflated
    water.expanded_volume = false; // mark the water as deflated

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = water.loc.x(), y = water.loc.y(), z = water.loc.z();
    int r = rh; // to avoid unsigned underflow
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) axes.x.bins); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) axes.y.bins); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) axes.z.bins); // zminus and zplus

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                // determine if the bin is within a sphere centered on the atom
                if (std::sqrt(std::pow(x - i, 2) + std::pow(y - j, 2) + std::pow(z - k, 2)) <= rh) {
                    if (grid.index(i, j, k) != GridObj::H_AREA) {continue;}
                    grid.index(i, j, k) = GridObj::EMPTY;
                }
            }
        }
    }
}

void Grid::clear_waters() {
    std::vector<Water> waters;
    waters.reserve(w_members.size());
    std::for_each(w_members.begin(), w_members.end(), [&waters] (const GridMember<Water>& water) {waters.push_back(water.atom);});
    remove(waters);
    if (w_members.size() != 0) [[unlikely]] {throw except::unexpected("Grid::clear_waters: Something went wrong.");}
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
    return std::pow(width, 3)*volume;
}

Grid& Grid::operator=(const Grid& rhs) {
    grid = rhs.grid;
    a_members = rhs.a_members;
    w_members = rhs.w_members;
    width = rhs.width;
    volume = rhs.volume;
    ra = rhs.ra;
    rh = rhs.rh;
    axes = rhs.axes;
    // culler & placer cannot be modified after program is run, so they'll automatically be equal always
    return *this;
}

Grid& Grid::operator=(Grid&& rhs) noexcept {
    grid = std::move(rhs.grid);
    a_members = std::move(rhs.a_members);
    w_members = std::move(rhs.w_members);
    width = rhs.width;
    volume = rhs.volume;
    ra = rhs.ra;
    rh = rhs.rh;
    axes = std::move(rhs.axes);
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

void Grid::save(const io::File& path) const {
    std::vector<Atom> atoms;
    std::vector<Water> waters;
    unsigned int c = 0;
    for (unsigned int i = 0; i < grid.xdim; i++) {
        for (unsigned int j = 0; j < grid.ydim; j++) {
            for (unsigned int k = 0; k < grid.zdim; k++) {
                switch (grid.index(i, j, k)) {
                    case GridObj::A_CENTER:
                        atoms.push_back(Atom(c++, "C", "", "LYS", "A", 1, "", Vector3<double>(i, j, k), 1, 0, constants::atom_t::C, ""));
                        break;
                    case GridObj::A_AREA:
                        atoms.push_back(Atom(c++, "C", "", "LYS", "B", 2, "", Vector3<double>(i, j, k), 1, 0, constants::atom_t::C, ""));
                        break;
                    case GridObj::H_CENTER:
                        waters.push_back(Water(c++, "H", "", "HOH", "C", 3, "", Vector3<double>(i, j, k), 1, 0, constants::atom_t::H, ""));
                        break;
                    case GridObj::H_AREA:
                        waters.push_back(Water(c++, "H", "", "HOH", "D", 4, "", Vector3<double>(i, j, k), 1, 0, constants::atom_t::H, ""));
                        break;
                }
            }
        }
    }
    Protein p(atoms, waters);
    p.save(path);
}

Body Grid::generate_excluded_volume() const {
    Body body;
    for (unsigned int i = 0; i < grid.xdim; i++) {
        for (unsigned int j = 0; j < grid.ydim; j++) {
            for (unsigned int k = 0; k < grid.zdim; k++) {
                switch (grid.index(i, j, k)) {
                    case GridObj::A_AREA:
                    case GridObj::A_CENTER: 
                        body.get_atoms().push_back(Water::create_new_water(to_xyz(i, j, k)));
                        break;
                    default: 
                        break;
                }
            }
        }
    }
    return body;
}

const GridObj::DATATYPE& Grid::index(unsigned int i, unsigned int j, unsigned int k) const {
    return grid.index(i, j, k);
}

std::vector<Atom> Grid::get_surface_atoms() const {
    throw except::not_implemented("Grid::get_surface_atoms: Not implemented.");
}

Vector3<double> Grid::to_xyz(int i, int j, int k) const {
    double x = axes.x.min + width*i;
    double y = axes.y.min + width*j;
    double z = axes.z.min + width*k;
    return {x, y, z};
}