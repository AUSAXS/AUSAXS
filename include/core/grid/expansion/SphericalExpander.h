// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include "grid/detail/GridObj.h"
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <utility/observer_ptr.h>
#include <settings/GridSettings.h>

namespace ausaxs::grid::volume {
    template<bool AtomicMinVol, bool WaterMinVol>
    struct SphericalExpander {
        static void expand_volume(observer_ptr<grid::Grid> grid, GridMember<data::AtomFF>& atom);
        static void expand_volume(observer_ptr<grid::Grid> grid, GridMember<data::Water>& atom);

        static void deflate_volume(observer_ptr<grid::Grid> grid, GridMember<data::AtomFF>& atom);
        static void deflate_volume(observer_ptr<grid::Grid> grid, GridMember<data::Water>& atom);
    };

    using AtomicExpander         = SphericalExpander<true, false>;
    using AtomicAndWaterExpander = SphericalExpander<true, true>;
}

template<bool AMV, bool _>
void ausaxs::grid::volume::SphericalExpander<AMV, _>::expand_volume(observer_ptr<grid::Grid> grid, GridMember<data::AtomFF>& atom) {
    if (atom.is_expanded()) {return;} // check if this location has already been expanded
    atom.set_expanded(true); // mark this location as expanded

    double rvdw = grid->get_atomic_radius(atom.get_atom_type());
    double rvdw2 = std::pow(rvdw, 2);
    [[maybe_unused]] double rvol = settings::grid::min_exv_radius;
    [[maybe_unused]] double rvol2 = std::pow(rvol, 2);

    int xm, xp, ym, yp, zm, zp;
    {   // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        auto axes = grid->get_axes();

        // determine maximum bin radius to check
        double br;
        if constexpr (AMV) {
            br = std::ceil(std::max(rvdw, rvol)/settings::grid::cell_width);
        } else {
            br = std::ceil(rvdw/settings::grid::cell_width);
        }

        auto [bx, by, bz] = atom.get_bin_loc();
        xm = std::max<int>(bx - br, 0), xp = std::min<int>(bx + br + 1, axes.x.bins); // xminus and xplus
        ym = std::max<int>(by - br, 0), yp = std::min<int>(by + br + 1, axes.y.bins); // yminus and yplus
        zm = std::max<int>(bz - br, 0), zp = std::min<int>(bz + br + 1, axes.z.bins); // zminus and zplus    
    }

    // loop over each bin in the box
    int added_volume = 0;

    // i, j, k *must* be ints to avoid unsigned underflow
    auto[ax, ay, az] = atom.get_absolute_loc();
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(grid->to_x(i) - ax, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(grid->to_y(j) - ay, 2);
            for (int k = zm; k < zp; ++k) {
                // fill a sphere of radius [0, vdw] around the atom
                double dist = x2y2 + std::pow(grid->to_z(k) - az, 2);
                auto& bin = grid->grid.index(i, j, k);

                bool fill = dist <= rvdw2;
                if (fill) {
                    if (!grid->grid.is_empty_or_volume(bin)) {continue;}
                    added_volume += !grid->grid.is_volume(bin); // only add to the volume if the bin is not already part of the volume
                    bin |= detail::A_AREA;
                }

                if constexpr (AMV) {
                    // fill outer shell of radius [vdw, rvol]
                    if (!fill && dist <= rvol2) {
                        if (!grid->grid.is_empty(i, j, k)) {continue;}
                        bin |= detail::VOLUME;
                        added_volume++;
                    }
                }
            }
        }
    }
    grid->add_volume(added_volume);
}

template<bool _, bool WMV>
void ausaxs::grid::volume::SphericalExpander<_, WMV>::expand_volume(observer_ptr<grid::Grid> grid, GridMember<data::Water>& water) {
    if (water.is_expanded()) {return;} // check if this location has already been expanded
    water.set_expanded(true); // mark this location as expanded

    double rvdw = grid->get_hydration_radius();
    double rvdw2 = std::pow(rvdw, 2);
    [[maybe_unused]] double rvol = settings::grid::min_exv_radius;
    [[maybe_unused]] double rvol2 = std::pow(rvol, 2);

    int xm, xp, ym, yp, zm, zp;
    {   // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        auto axes = grid->get_axes();

        // determine maximum bin radius to check
        double br;
        if constexpr (WMV) {
            br = std::ceil(std::max(rvdw, rvol)/settings::grid::cell_width);
        } else {
            br = std::ceil(rvdw/settings::grid::cell_width);
        }

        auto[bx, by, bz] = water.get_bin_loc();
        xm = std::max<int>(bx - br, 0), xp = std::min<int>(bx + br + 1, axes.x.bins); // xminus and xplus
        ym = std::max<int>(by - br, 0), yp = std::min<int>(by + br + 1, axes.y.bins); // yminus and yplus
        zm = std::max<int>(bz - br, 0), zp = std::min<int>(bz + br + 1, axes.z.bins); // zminus and zplus    
    }

    // i, j, k *must* be ints to avoid unsigned underflow
    auto [ax, ay, az] = water.get_absolute_loc();
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(grid->to_x(i) - ax, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(grid->to_y(j) - ay, 2);
            for (int k = zm; k < zp; ++k) {
                // determine if the bin is within a sphere centered on the atom
                double dist = x2y2 + std::pow(grid->to_z(k) - az, 2);
                
                bool fill = dist <= rvdw2;
                if (fill) {
                    if (!grid->grid.is_empty(i, j, k)) {continue;}
                    grid->grid.index(i, j, k) |= detail::W_AREA;
                }

                if constexpr (WMV) {
                    // fill outer shell of radius [vdw, rvol]
                    if (!fill && dist <= rvol2) {
                        if (!grid->grid.is_empty(i, j, k)) {continue;}
                        grid->grid.index(i, j, k) |= detail::VOLUME;
                    }
                }
            }
        }
    }
}

template<bool AMV, bool _>
void ausaxs::grid::volume::SphericalExpander<AMV, _>::deflate_volume(observer_ptr<grid::Grid> grid, GridMember<data::AtomFF>& atom) {
    if (!atom.is_expanded()) {return;} // check if this location has already been deflated
    atom.set_expanded(false); // mark the atom as deflated

    double rmax;
    if constexpr (AMV) {
        rmax = std::max(grid->get_atomic_radius(atom.get_atom_type()), settings::grid::min_exv_radius);
    } else {
        rmax = grid->get_atomic_radius(atom.get_atom_type());
    }
    double rmax2 = std::pow(rmax, 2);

    int xm, xp, ym, yp, zm, zp;
    {   // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        auto axes = grid->get_axes();
        double br = std::ceil(rmax/settings::grid::cell_width);
        auto [bx, by, bz] = atom.get_bin_loc();
        xm = std::max<int>(bx - br, 0), xp = std::min<int>(bx + br + 1, axes.x.bins); // xminus and xplus
        ym = std::max<int>(by - br, 0), yp = std::min<int>(by + br + 1, axes.y.bins); // yminus and yplus
        zm = std::max<int>(bz - br, 0), zp = std::min<int>(bz + br + 1, axes.z.bins); // zminus and zplus    
    }

    // i, j, k *must* be ints due to avoid unsigned underflow
    int removed_volume = 0;
    auto [ax, ay, az] = atom.get_absolute_loc();
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(grid->to_x(i) - ax, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(grid->to_y(j) - ay, 2);
            for (int k = zm; k < zp; ++k) {
                double dist = x2y2 + std::pow(grid->to_z(k) - az, 2);

                // determine if the bin is within a sphere centered on the atom
                auto& bin = grid->grid.index(i, j, k);
                if (dist <= rmax2) {
                    removed_volume += grid->grid.is_only_atom_area_or_volume(bin);
                    bin &= ~(detail::A_AREA | detail::VOLUME);
                }
            }
        }
    }
    grid->add_volume(-removed_volume); // only the actual atoms contributes to the volume
}

template<bool _, bool WMV>
void ausaxs::grid::volume::SphericalExpander<_, WMV>::deflate_volume(observer_ptr<grid::Grid> grid, GridMember<data::Water>& water) {
    if (!water.is_expanded()) {return;} // check if this location has already been deflated
    water.set_expanded(false); // mark the water as deflated

    double rmax;
    if constexpr (WMV) {
        rmax = std::max(grid->get_hydration_radius(), settings::grid::min_exv_radius);
    } else {
        rmax = grid->get_hydration_radius();
    }
    double rmax2 = std::pow(rmax, 2);

    int xm, xp, ym, yp, zm, zp;
    {   // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
        auto axes = grid->get_axes();
        auto[bx, by, bz] = water.get_bin_loc();
        double r = std::ceil(rmax/settings::grid::cell_width);
        xm = std::max<int>(bx - r, 0), xp = std::min<int>(bx + r + 1, axes.x.bins); // xminus and xplus
        ym = std::max<int>(by - r, 0), yp = std::min<int>(by + r + 1, axes.y.bins); // yminus and yplus
        zm = std::max<int>(bz - r, 0), zp = std::min<int>(bz + r + 1, axes.z.bins); // zminus and zplus    
    }

    // i, j, k *must* be ints to avoid unsigned underflow
    auto [ax, ay, az] = water.get_absolute_loc();
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(grid->to_x(i) - ax, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(grid->to_y(j) - ay, 2);
            for (int k = zm; k < zp; ++k) {
                double dist = x2y2 + std::pow(grid->to_z(k) - az, 2);
                if (dist <= rmax2) {
                    grid->grid.index(i, j, k) &= ~detail::W_AREA;
                    // grid->grid.index(i, j, k) &= ~(detail::W_AREA | detail::VOLUME);
                }
            }
        }
    }
}