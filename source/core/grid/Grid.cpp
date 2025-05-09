/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "grid/detail/GridObj.h"
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <grid/detail/GridSurfaceDetection.h>
#include <grid/expansion/GridExpander.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/atoms/AtomHelper.h>
#include <settings/GridSettings.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <constants/Constants.h>
#include <io/ExistingFile.h>

using namespace ausaxs;
using namespace ausaxs::grid;
using namespace ausaxs::data;

Grid::Grid(const Limit3D& axes) : axes(Axis3D(axes, settings::grid::cell_width)) {
    setup();
}

Grid::Grid(const std::vector<AtomFF>& atoms) : Grid({Body(atoms)}) {}

Grid::Grid(const std::vector<Body>& bodies) {
    // find the total bounding box containing all bodies
    Vector3 min{0, 0, 0}, max{0, 0, 0};
    for (const Body& body : bodies) {
        auto[amin, amax] = bounding_box(body.get_atoms());

        for (unsigned int i = 0; i < 3; i++) {
            if (amin[i] < min[i]) {min[i] = amin[i];}
            if (amax[i] > max[i]) {max[i] = amax[i];}
        }

        auto w = body.get_waters();
        if (!w.has_value()) {continue;}
        auto[wmin, wmax] = bounding_box(w.value().get());
        for (unsigned int i = 0; i < 3; i++) {
            if (wmin[i] < min[i]) {min[i] = wmin[i];}
            if (wmax[i] > max[i]) {max[i] = wmax[i];}
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
        nmin[i] = std::floor(min[i] - expand - settings::grid::cell_width); // flooring to make our grid 'nice' (i.e. bin edges are at integer values)
        nmax[i] = std::ceil( max[i] + expand + settings::grid::cell_width); //  ceiling to make our grid 'nice' (i.e. bin edges are at integer values)
    }

    // setup the rest of the class members
    axes = Axis3D(nmin, nmax, settings::grid::cell_width);
    setup();

    // finally add all atoms to the grid
    for (const Body& body : bodies) {
        add(body, false);
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
        double min_length = 0.5*settings::grid::min_bins*settings::grid::cell_width;
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

double Grid::get_atomic_radius(form_factor::form_factor_t atom) const {
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
            if (min[i] > atom.get_bin_loc()[i]) min[i] = atom.get_bin_loc()[i];     // min
            if (max[i] < atom.get_bin_loc()[i]) max[i] = atom.get_bin_loc()[i]+1;   // max. +1 since this will often be used as loop limits
        }
    }
    return std::make_pair(min, max);
}

template<data::AtomType T>
std::pair<Vector3<double>, Vector3<double>> _bounding_box(const std::vector<T>& atoms) {
    // initialize the bounds as extreme as possible
    Vector3 min = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vector3 max = {std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min()};
    for (const auto& atom : atoms) {
        for (int i = 0; i < 3; i++) {
            min[i] = std::min(min[i], atom.coordinates()[i]);
            max[i] = std::max(max[i], atom.coordinates()[i]);
        }
    }
    return std::make_pair(min, max);
}

std::pair<Vector3<double>, Vector3<double>> Grid::bounding_box(const std::vector<AtomFF>& atoms) {
    return _bounding_box<AtomFF>(atoms);
}

std::pair<Vector3<double>, Vector3<double>> Grid::bounding_box(const std::vector<Water>& atoms) {
    return _bounding_box<Water>(atoms);
}

void Grid::add_volume(double value) {
    volume += value;
}

void Grid::force_expand_volume() {
    for (auto& atom : a_members) {
        atom.set_expanded(false);
        volume::expand(this, atom);
    }

    for (auto& water : w_members) {
        water.set_expanded(false);
        volume::expand(this, water);
    }
}

void Grid::expand_volume() {
    // iterate through each member location
    for (auto& atom : a_members) {
        volume::expand(this, atom);
    }

    for (auto& water : w_members) {
        volume::expand(this, water);
    }
}

void Grid::deflate_volume() {
    // iterate through each member location
    for (auto& atom : a_members) {
        volume::deflate(this, atom);
    }

    for (auto& water : w_members) {
        volume::deflate(this, water);
    }
}

void Grid::remove_waters(const std::vector<bool>& to_remove) {
    assert(
        to_remove.size() == w_members.size() && 
        "Grid::remove_waters: The size of the removal vector does not match the number of waters!"
    );

    std::vector<GridMember<Water>> new_waters;
    new_waters.reserve(w_members.size());
    for (unsigned int i = 0; i < w_members.size(); ++i) {
        if (!to_remove[i]) {
            new_waters.push_back(w_members[i]);
        } else {
            volume::deflate(this, w_members[i]);

            // clear only W_CENTER since bin may be shared with A_CENTER
            grid.index(w_members[i].get_bin_loc()) = detail::EMPTY;
        }
    }
    w_members = std::move(new_waters);
}

std::span<GridMember<AtomFF>> Grid::add(const Body& body, bool expand) {
    int start = a_members.size();
    body_start[body.get_uid()] = start;
    if (body.size_atom() == 0) {
        return {a_members.begin(), a_members.end()};
    }

    a_members.resize(a_members.size() + body.symmetry().size_atom_total());
    auto b_atoms = body.symmetry().explicit_structure().atoms;

    for (int i = start; i < static_cast<int>(a_members.size()); i++) {
        auto& atom = b_atoms[i-start];
        auto loc = to_bins(atom.coordinates());
        unsigned int x = loc.x(), y = loc.y(), z = loc.z();

        // sanity check
        #if DEBUG
            bool out_of_bounds = x >= axes.x.bins || y >= axes.y.bins || z >= axes.z.bins;
            if (out_of_bounds) [[unlikely]] {
                throw except::out_of_bounds(
                    "Grid::add: Atom is located outside the grid!\nBin location: "
                     + loc.to_string() + "\n: " + axes.to_string() + "\n"
                    "Real location: " + atom.coordinates().to_string()
                );
            }
        #endif

        auto& bin = grid.index(x, y, z);
        volume += grid.is_empty_or_water(bin);
        bin = detail::A_CENTER;

        GridMember gm(atom, std::move(loc));
        if (expand) {volume::expand(this, gm);}
        a_members[i] = std::move(gm);
    }

    // add the hydration if present
    auto w = body.get_waters();
    if (w.has_value()) {
        add(w.value().get(), expand);
    }

    return {a_members.begin() + start, a_members.end()};
}

auto add_single_water = [] (grid::Grid& g, const data::Water& w) {
    auto loc = g.to_bins(w.coordinates());
    int x = loc.x(), y = loc.y(), z = loc.z(); 

    // sanity check
    #if DEBUG
        auto ax = g.get_axes();
        if (x >= (int) ax.x.bins || y >= (int) ax.y.bins || z >= (int) ax.z.bins || x < 0 || y < 0 || z < 0) {
            throw except::out_of_bounds(
                "Grid::add: Atom is located outside the grid!\nBin location: " + loc.to_string() + "\n: " + g.get_axes().to_string() + "\n"
                "Real location: " + w.coordinates().to_string()
            );
        }

        auto bin = g.index(x, y, z);
        if (!(bin == grid::detail::EMPTY || bin == grid::detail::W_AREA || bin == grid::detail::VOLUME)) {
            throw except::invalid_operation(
                "Grid::add: Attempting to add a water molecule to a non-empty location containing " + std::to_string(bin) +
                " (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")"
            );
        }
    #endif
    g.grid.index(x, y, z) = grid::detail::W_CENTER*g.grid.is_empty_or_volume(x, y, z);
    return GridMember(w, std::move(loc));
};

std::span<grid::GridMember<data::Water>> Grid::add(const std::vector<data::Water>& waters, bool expand) {
    size_t start = w_members.size();
    w_members.reserve(w_members.size() + waters.size());
    for (int i = 0; i < static_cast<int>(waters.size()); ++i) {
        auto& w = waters[i];
        auto gm = add_single_water(*this, w);
        if (expand) {volume::expand(this, gm);}
        w_members.emplace_back(std::move(gm));
    }
    return {w_members.begin() + start, w_members.end()};
}

grid::GridMember<data::Water>& Grid::add(const data::Water& water, bool expand) {
    w_members.emplace_back(add_single_water(*this, water));
    if (expand) {volume::expand(this, w_members.back());}
    return w_members.back();
}

void Grid::remove(const Body& body) {
    assert(body_start.contains(body.get_uid()) && "Grid::remove: Attempting to remove a body that is not in the grid!");
    auto it = body_start.find(body.get_uid());

    // find start and end indices of the body atoms
    int start = it->second;
    int end = it->second + body.symmetry().size_atom_total();
    assert(end <= static_cast<int>(a_members.size()) && "Grid::remove: Contained bodies have been modified after being added to the grid.");

    // deflate and remove all atoms in the body from the grid
    for (int i = start; i < end; i++) {
        volume::deflate(this, a_members[i]);
        auto& bin = grid.index(a_members[i].get_bin_loc());
        volume -= !grid.is_empty(bin); // multiple atoms may share a center bin, so we have to check if its volume was already removed
        bin = detail::EMPTY;
    }

    int diff = end - start;

    // clean up the internal data
    a_members.erase(a_members.begin()+start, a_members.begin()+end); // erase its atoms
    body_start.erase(it);                   // erase removed body
    for (auto& [_, start] : body_start) {   // update start indices of remaining bodies
        if (end <= start) {
            start -= diff;
        }
    }
}

void Grid::clear_waters() {
    for (auto& water : w_members) {
        volume::deflate(this, water);
        grid.index(water.get_bin_loc()) = detail::EMPTY;
    }
    w_members.clear();
}

Vector3<int> Grid::get_bins() const {
    return Vector3<int>(axes.x.bins, axes.y.bins, axes.z.bins);
}

Vector3<int> Grid::to_bins(const Vector3<double>& v) const {
    int binx = std::round((v.x() - axes.x.min)/settings::grid::cell_width);
    int biny = std::round((v.y() - axes.y.min)/settings::grid::cell_width);
    int binz = std::round((v.z() - axes.z.min)/settings::grid::cell_width);
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
    return volume*std::pow(settings::grid::cell_width, 3);
}

Grid& Grid::operator=(const Grid& rhs) {
    grid = rhs.grid;
    a_members = rhs.a_members;
    w_members = rhs.w_members;
    volume = rhs.volume;
    axes = rhs.axes;
    body_start = rhs.body_start;
    return *this;
}

Grid& Grid::operator=(Grid&& rhs) noexcept {
    grid = std::move(rhs.grid);
    a_members = std::move(rhs.a_members);
    w_members = std::move(rhs.w_members);
    volume = rhs.volume;
    axes = std::move(rhs.axes);
    body_start = std::move(rhs.body_start);
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
    std::vector<std::vector<AtomFF>> atoms(7);
    for (unsigned int i = 0; i < grid.size_x(); i++) {
        for (unsigned int j = 0; j < grid.size_y(); j++) {
            for (unsigned int k = 0; k < grid.size_z(); k++) {
                switch (grid.index(i, j, k)) {
                    case detail::A_CENTER:
                        atoms[0].emplace_back(to_xyz(i, j, k), form_factor::form_factor_t::C);
                        break;
                    case detail::A_AREA:
                        atoms[1].emplace_back(to_xyz(i, j, k), form_factor::form_factor_t::C);
                        break;
                    case detail::VOLUME:
                        atoms[2].emplace_back(to_xyz(i, j, k), form_factor::form_factor_t::C);
                        break;
                    case detail::W_CENTER:
                        atoms[3].emplace_back(to_xyz(i, j, k), form_factor::form_factor_t::C);
                        break;
                    case detail::W_AREA:
                        atoms[4].emplace_back(to_xyz(i, j, k), form_factor::form_factor_t::C);
                        break;
                    case detail::VACUUM:
                        atoms[5].emplace_back(to_xyz(i, j, k), form_factor::form_factor_t::C);
                    default:
                        break;
                }
            }
        }
    }

    // visualize corners
    atoms[6].emplace_back(to_xyz(0,             0,              0),           form_factor::form_factor_t::C);
    atoms[6].emplace_back(to_xyz(axes.x.bins,   0,              0),           form_factor::form_factor_t::C);
    atoms[6].emplace_back(to_xyz(axes.x.bins,   axes.y.bins,    0),           form_factor::form_factor_t::C);
    atoms[6].emplace_back(to_xyz(0,             axes.y.bins,    0),           form_factor::form_factor_t::C);
    atoms[6].emplace_back(to_xyz(0,             axes.y.bins,    axes.z.bins), form_factor::form_factor_t::C);
    atoms[6].emplace_back(to_xyz(0,             0,              axes.z.bins), form_factor::form_factor_t::C);
    atoms[6].emplace_back(to_xyz(axes.x.bins,   0,              axes.z.bins), form_factor::form_factor_t::C);
    atoms[6].emplace_back(to_xyz(axes.x.bins,   axes.y.bins,    axes.z.bins), form_factor::form_factor_t::C);
    
    std::vector<Body> bodies;
    bodies.emplace_back(atoms[0]);
    bodies.emplace_back(atoms[1]);
    bodies.emplace_back(atoms[2]);
    bodies.emplace_back(atoms[3]);
    bodies.emplace_back(atoms[4]);
    bodies.emplace_back(atoms[5]);
    bodies.emplace_back(atoms[6]);

    data::Molecule p(bodies);
    p.save(path);
}

grid::detail::GridExcludedVolume Grid::generate_excluded_volume(bool determine_surface) {
    expand_volume();
    auto vol = determine_surface ? detail::GridSurfaceDetection(this).detect() : detail::GridSurfaceDetection(this).no_detect();

    if (settings::grid::exv::save) {
        std::vector<AtomFF> atoms1, atoms2;
    
        for (int i = 0; i < static_cast<int>(vol.interior.size()); ++i) {
            atoms1.emplace_back(vol.interior[i], form_factor::form_factor_t::C);
        }

        for (int j = 0; j < static_cast<int>(vol.surface.size()); ++j) {
            atoms2.emplace_back(vol.surface[j], form_factor::form_factor_t::C);
        }

        std::vector<Body> bodies = {Body{atoms1}, Body{atoms2}};
        data::Molecule(bodies).save(settings::general::output + "exv.pdb");
    }

    if (!determine_surface) {
        return {vol.interior, {}};
    }
    return vol;
}

std::vector<data::AtomFF> Grid::get_atoms() {
    std::vector<data::AtomFF> atoms;
    for (auto& atom : a_members) {
        atoms.emplace_back(atom.get_atom());
    }
    return atoms;
}

std::vector<data::Water> Grid::get_waters() {
    std::vector<data::Water> waters;
    for (auto& water : w_members) {
        waters.emplace_back(water.get_atom());
    }
    return waters;
}

const grid::detail::State& Grid::index(unsigned int i, unsigned int j, unsigned int k) const {
    return grid.index(i, j, k);
}

std::vector<AtomFF> Grid::get_surface_atoms() const {
    throw except::not_implemented("Grid::get_surface_atoms: Not implemented.");
}

Vector3<double> Grid::to_xyz(int i, int j, int k) const {
    return {to_x(i), to_y(j), to_z(k)};
}

double Grid::to_x(int i) const {
    return axes.x.min + i*settings::grid::cell_width;
}

double Grid::to_y(int j) const {
    return axes.y.min + j*settings::grid::cell_width;
}

double Grid::to_z(int k) const {
    return axes.z.min + k*settings::grid::cell_width;
}

Vector3<int> Grid::get_center() const {
    return {int(axes.x.bins/2), int(axes.y.bins/2), int(axes.z.bins/2)};
}

double Grid::get_width() const {return settings::grid::cell_width;}