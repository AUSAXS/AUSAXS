/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/detail/GridSurfaceDetection.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/record/Atom.h>
#include <settings/GridSettings.h>

using namespace grid::detail;

GridSurfaceDetection::GridSurfaceDetection(observer_ptr<grid::Grid> grid) : RadialLineGenerator(grid, {std::sqrt(grid->get_width())+1e-3, 2*grid->get_width(), 3*grid->get_width(), 4*grid->get_width()}), grid(grid) {}

GridSurfaceDetection::~GridSurfaceDetection() = default;

bool GridSurfaceDetection::collision_check(const Vector3<int>& loc) const {
    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    int score = 0;

    // check if a location is out-of-bounds
    auto is_out_of_bounds = [&bins] (Vector3<int> v) {
        if (v.x() < 0 || (int) bins.x() <= v.x()) {return true;}
        if (v.y() < 0 || (int) bins.y() <= v.y()) {return true;}
        if (v.z() < 0 || (int) bins.z() <= v.z()) {return true;}
        return false;
    };

    for (unsigned int i = 0; i < rot_locs_abs.size(); i++) {
        {   // check for collisions at 1r
            auto xr = loc.x() + rot_bins_1[i].x();
            auto yr = loc.y() + rot_bins_1[i].y();
            auto zr = loc.z() + rot_bins_1[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 7;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                continue;
            }
        }

        {   // check for collisions at 2r
            auto xr = loc.x() + rot_bins_2[i].x();
            auto yr = loc.y() + rot_bins_2[i].y();
            auto zr = loc.z() + rot_bins_2[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 7;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                score += 3;
                continue;
            }
        }

        {   // check for collisions at 3r
            auto xr = loc.x() + rot_bins_3[i].x();
            auto yr = loc.y() + rot_bins_3[i].y();
            auto zr = loc.z() + rot_bins_3[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 7;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                score += 5;
                continue;
            }
        }

        score += 7;
    }
    return score < 42;
}

bool GridSurfaceDetection::vacuum_collision_check(const Vector3<int>& loc) const {
    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    // check if a location is out-of-bounds
    auto is_out_of_bounds = [&bins] (Vector3<int> v) {
        if (v.x() < 0 || (int) bins.x() <= v.x()) {return true;}
        if (v.y() < 0 || (int) bins.y() <= v.y()) {return true;}
        if (v.z() < 0 || (int) bins.z() <= v.z()) {return true;}
        return false;
    };

    int count = 0;
    for (unsigned int i = 0; i < rot_locs_abs.size(); ++i) {
        // check for collisions at 1r
        auto xr = loc.x() + rot_bins_1[i].x();
        auto yr = loc.y() + rot_bins_1[i].y();
        auto zr = loc.z() + rot_bins_1[i].z();
        if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
            return false; // if we are out of bounds we must be on the surface
        }

        if (!gref.is_empty_or_water(xr, yr, zr)) {
            ++count; // collision at 1r means the hole must be too small for a water molecule
        }
    }
    return count <= 1; // allow op to one collision at 1r with the parent voxel
}

std::vector<Vector3<double>> GridSurfaceDetection::determine_vacuum_holes(const std::vector<Vector3<int>>& surface) const {
    std::vector<Vector3<double>> vacuum_voxels; vacuum_voxels.reserve(surface.size());
    for (const auto& loc : surface) {
        for (unsigned int i = 0; i < rot_locs_abs.size(); ++i) {
            auto xr = loc.x() + rot_bins_1[i].x();
            auto yr = loc.y() + rot_bins_1[i].y();
            auto zr = loc.z() + rot_bins_1[i].z();
            if (grid->grid.is_empty_or_water(xr, yr, zr) && vacuum_collision_check({xr, yr, zr})) {
                vacuum_voxels.push_back(grid->to_xyz(xr, yr, zr));
            }
        }
    }

    // fill free space around vacuum voxels
    int stride = std::round(2*settings::grid::exv_radius/settings::grid::width);
    int buffer = 1./settings::grid::width; // 2Å cube should be enough to be space-filling
    const auto& axes = grid->get_axes();
    auto& gobj = grid->grid;
    for (const auto& loc : vacuum_voxels) {
        for (int i = std::max<int>(loc.x()-buffer, 0); i < std::min<int>(loc.x()+buffer, axes.x.bins); i+=stride) {
            for (int j = std::max<int>(loc.y()-buffer, 0); j < std::min<int>(loc.y()+buffer, axes.y.bins); j+=stride) {
                for (int k = std::max<int>(loc.z()-buffer, 0); k < std::min<int>(loc.z()+buffer, axes.z.bins); k+=stride) {
                    if (gobj.is_empty_or_water(i, j, k)) {
                        gobj.index(i, j, k) = grid::detail::VACUUM;
                        vacuum_voxels.push_back(grid->to_xyz(i, j, k));
                    }
                }
            }
        }
    }

    return vacuum_voxels;
}

template<bool detect_surface>
GridExcludedVolume GridSurfaceDetection::helper() const {
    GridExcludedVolume vol;
    vol.interior.reserve(grid->get_volume());

    int stride = std::round(2*settings::grid::exv_radius/settings::grid::width);
    int buffer = std::round(std::max<double>(settings::grid::rvol, 2)/settings::grid::width);
    const auto& axes = grid->get_axes();
    const auto& gobj = grid->grid;
    auto[vmin, vmax] = grid->bounding_box_index();
    for (int i = std::max<int>(vmin.x()-buffer, 0); i < std::min<int>(vmax.x()+buffer+1, axes.x.bins); i+=stride) {
        for (int j = std::max<int>(vmin.y()-buffer, 0); j < std::min<int>(vmax.y()+buffer+1, axes.y.bins); j+=stride) {
            for (int k = std::max<int>(vmin.z()-buffer, 0); k < std::min<int>(vmax.z()+buffer+1, axes.z.bins); k+=stride) {
                switch (gobj.index(i, j, k)) {
                    case grid::detail::VOLUME:
                    case grid::detail::A_AREA:
                    case grid::detail::A_CENTER: 
                        if constexpr (detect_surface) {
                            if (collision_check({i, j, k})) {
                                vol.interior.push_back(grid->to_xyz(i, j, k));
                            } else {
                                vol.surface.push_back(grid->to_xyz(i, j, k));
                            }
                        } else {
                            vol.interior.push_back(grid->to_xyz(i, j, k));
                        }
                    default:
                        break;
                }
            }
        }
    }

    std::vector<Vector3<int>> surface;
    std::transform(vol.surface.begin(), vol.surface.end(), std::back_inserter(surface), [this] (const Vector3<double>& v) {return grid->to_bins(v);});
    vol.vacuum = determine_vacuum_holes(surface);

    return vol;
}
template GridExcludedVolume GridSurfaceDetection::helper<true>() const;
template GridExcludedVolume GridSurfaceDetection::helper<false>() const;

GridExcludedVolume GridSurfaceDetection::no_detect() const {
    return helper<false>();
}

GridExcludedVolume GridSurfaceDetection::detect() const {
    return helper<true>();
}