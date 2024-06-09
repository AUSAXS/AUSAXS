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
                score += 3;
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
                score += 3;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                score += 1;
                continue;
            }
        }

        {   // check for collisions at 3r
            auto xr = loc.x() + rot_bins_3[i].x();
            auto yr = loc.y() + rot_bins_3[i].y();
            auto zr = loc.z() + rot_bins_3[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 3;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                score += 2;
                continue;
            }
        }

        score += 3;
    }
    return score < 14;
}

template<bool detect_surface>
GridExcludedVolume GridSurfaceDetection::helper() const {
    GridExcludedVolume vol;
    vol.interior.reserve(grid->get_volume());

    int stride = std::round(2*settings::grid::exv_radius/settings::grid::width);
    int buffer = 2./settings::grid::width; // 2Å buffer in each direction should be enough to capture all filled voxels
    const auto& axes = grid->get_axes();
    const auto& gobj = grid->grid;
    auto[vmin, vmax] = grid->bounding_box_index();
    for (int i = std::max<int>(vmin.x()-buffer, 0); i < std::min<int>(vmax.x()+buffer, axes.x.bins); i+=stride) {
        for (int j = std::max<int>(vmin.y()-buffer, 0); j < std::min<int>(vmax.y()+buffer, axes.y.bins); j+=stride) {
            for (int k = std::max<int>(vmin.z()-buffer, 0); k < std::min<int>(vmax.z()+buffer, axes.z.bins); k+=stride) {
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