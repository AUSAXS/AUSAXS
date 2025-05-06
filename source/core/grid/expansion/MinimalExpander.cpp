/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/expansion/MinimalExpander.h>
#include <grid/Grid.h>
#include <settings/GridSettings.h>

using namespace ausaxs;
using namespace ausaxs::grid::expander;

void MinimalExpander::expand_volume(GridMember<data::AtomFF>& atom) {
    if (atom.is_expanded()) {return;} // check if this location has already been expanded
    atom.set_expanded(true); // mark this location as expanded

    auto[ax, ay, az] = atom.get_absolute_loc();
    double rvdw2 = std::pow(grid->get_atomic_radius(atom.get_atom_type()), 2);
    double rvol2 = std::pow(settings::grid::min_exv_radius, 2);

    int xm, xp, ym, yp, zm, zp;
    {   // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        auto axes = grid->get_axes();
        double br = std::ceil(std::max(grid->get_atomic_radius(atom.get_atom_type()), settings::grid::min_exv_radius)/settings::grid::cell_width);
        auto [bx, by, bz] = atom.get_bin_loc();
        xm = std::max<int>(bx - br, 0), xp = std::min<int>(bx + br + 1, axes.x.bins); // xminus and xplus
        ym = std::max<int>(by - br, 0), yp = std::min<int>(by + br + 1, axes.y.bins); // yminus and yplus
        zm = std::max<int>(bz - br, 0), zp = std::min<int>(bz + br + 1, axes.z.bins); // zminus and zplus    
    }

    // loop over each bin in the box
    int added_volume = 0;

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(to_x(i) - ax, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(to_y(j) - ay, 2);
            for (int k = zm; k < zp; ++k) {
                // fill a sphere of radius [0, vdw] around the atom
                double dist = x2y2 + std::pow(to_z(k) - az, 2);
                auto& bin = grid->grid.index(i, j, k);
                if (dist <= rvdw2) {
                    if (!grid->grid.is_empty_or_volume(bin)) {continue;}
                    added_volume += !grid->grid.is_volume(bin); // only add to the volume if the bin is not already part of the volume
                    bin = detail::A_AREA;
                }

                // fill an outer shell of radius [vdw, rvol] to make sure the volume is space-filling
                else if (dist <= rvol2) {
                    if (!grid->grid.is_empty(i, j, k)) {continue;}
                    bin = detail::VOLUME;
                    added_volume++;
                }
            }
        }
    }
    get_volume() += added_volume;
}

void MinimalExpander::expand_volume(GridMember<data::Water>& water) {
    if (water.is_expanded()) {return;} // check if this location has already been expanded
    water.set_expanded(true); // mark this location as expanded

    auto [ax, ay, az] = water.get_absolute_loc();
    double rvdw2 = std::pow(grid->get_hydration_radius(), 2);

    int xm, xp, ym, yp, zm, zp;
    {   // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        auto axes = grid->get_axes();
        auto[bx, by, bz] = water.get_bin_loc();
        double br = std::ceil(grid->get_hydration_radius()/settings::grid::cell_width);
        xm = std::max<int>(bx - br, 0), xp = std::min<int>(bx + br + 1, axes.x.bins); // xminus and xplus
        ym = std::max<int>(by - br, 0), yp = std::min<int>(by + br + 1, axes.y.bins); // yminus and yplus
        zm = std::max<int>(bz - br, 0), zp = std::min<int>(bz + br + 1, axes.z.bins); // zminus and zplus    
    }

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(to_x(i) - ax, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(to_y(j) - ay, 2);
            for (int k = zm; k < zp; ++k) {
                // determine if the bin is within a sphere centered on the atom
                double dist = x2y2 + std::pow(to_z(k) - az, 2);
                if (dist <= rvdw2) {
                    if (!grid->grid.is_empty(i, j, k)) {continue;}
                    grid->grid.index(i, j, k) = detail::W_AREA;
                }
            }
        }
    }
}

void MinimalExpander::deflate_volume(GridMember<data::AtomFF>& atom) {
    if (!atom.is_expanded()) {return;} // check if this location has already been deflated
    atom.set_expanded(false); // mark the atom as deflated

    auto [ax, ay, az] = atom.get_absolute_loc();
    double rvol = std::max(grid->get_atomic_radius(atom.get_atom_type()), settings::grid::min_exv_radius);
    double rvol2 = std::pow(rvol, 2);

    int xm, xp, ym, yp, zm, zp;
    {   // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        auto axes = grid->get_axes();
        auto [bx, by, bz] = atom.get_bin_loc();
        double br = std::ceil(rvol/settings::grid::cell_width);
        xm = std::max<int>(bx - br, 0), xp = std::min<int>(bx + br + 1, axes.x.bins); // xminus and xplus
        ym = std::max<int>(by - br, 0), yp = std::min<int>(by + br + 1, axes.y.bins); // yminus and yplus
        zm = std::max<int>(bz - br, 0), zp = std::min<int>(bz + br + 1, axes.z.bins); // zminus and zplus    
    }

    // i, j, k *must* be ints due to avoid unsigned underflow
    int removed_volume = 0;
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(to_x(i) - ax, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(to_y(j) - ay, 2);
            for (int k = zm; k < zp; ++k) {
                double dist = x2y2 + std::pow(to_z(k) - az, 2);

                // determine if the bin is within a sphere centered on the atom
                auto& bin = grid->grid.index(i, j, k);
                if (dist <= rvol2) {
                    if (!grid->grid.is_atom_area_or_volume(bin)) {continue;}
                    bin = detail::EMPTY;
                    ++removed_volume;
                }
            }
        }
    }
    get_volume() -= removed_volume; // only the actual atoms contributes to the volume
}

void MinimalExpander::deflate_volume(GridMember<data::Water>& water) {
    if (!water.is_expanded()) {return;} // check if this location has already been deflated
    water.set_expanded(false); // mark the water as deflated

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    auto [ax, ay, az] = water.get_absolute_loc();
    double rvdw = grid->get_hydration_radius();
    double rvdw2 = std::pow(rvdw, 2);

    int xm, xp, ym, yp, zm, zp;
    {   // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        auto axes = grid->get_axes();
        auto[bx, by, bz] = water.get_bin_loc();
        double r = std::ceil(rvdw/settings::grid::cell_width);
        xm = std::max<int>(bx - r, 0), xp = std::min<int>(bx + r + 1, axes.x.bins); // xminus and xplus
        ym = std::max<int>(by - r, 0), yp = std::min<int>(by + r + 1, axes.y.bins); // yminus and yplus
        zm = std::max<int>(bz - r, 0), zp = std::min<int>(bz + r + 1, axes.z.bins); // zminus and zplus    
    }

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(to_x(i) - ax, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(to_y(j) - ay, 2);
            for (int k = zm; k < zp; ++k) {
                double dist = x2y2 + std::pow(to_z(k) - az, 2);
                if (dist <= rvdw2) {
                    if (!grid->grid.is_water_area(i, j, k)) {continue;}
                    grid->grid.index(i, j, k) = detail::EMPTY;
                }
            }
        }
    }
}