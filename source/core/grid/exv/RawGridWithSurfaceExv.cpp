#include <grid/exv/RawGridWithSurfaceExv.h>
#include <grid/detail/RadialLineGenerator.h>
#include <grid/Grid.h>
#include <grid/detail/GridObj.h>
#include <utility/Logging.h>
#include <settings/GridSettings.h>

using namespace ausaxs;
using namespace ausaxs::grid::exv;

bool collision_check(observer_ptr<grid::Grid> grid, const Vector3<int>& loc) {
    static double width = grid->get_width();
    static grid::detail::RadialLineGenerator radial_lines(grid, {std::sqrt(grid->get_width())+1e-3, 2*grid->get_width(), 3*grid->get_width(), 4*grid->get_width()});
    if (width != grid->get_width()) {
        width = grid->get_width();
        radial_lines = grid::detail::RadialLineGenerator(grid, {std::sqrt(grid->get_width())+1e-3, 2*grid->get_width(), 3*grid->get_width(), 4*grid->get_width()});
    }

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

    for (int i = 0; i < static_cast<int>(radial_lines.rot_locs_abs.size()); i++) {
        {   // check for collisions at 1r
            auto xr = loc.x() + radial_lines.rot_bins_1[i].x();
            auto yr = loc.y() + radial_lines.rot_bins_1[i].y();
            auto zr = loc.z() + radial_lines.rot_bins_1[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 7;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                continue;
            }
        }

        {   // check for collisions at 2r
            auto xr = loc.x() + radial_lines.rot_bins_2[i].x();
            auto yr = loc.y() + radial_lines.rot_bins_2[i].y();
            auto zr = loc.z() + radial_lines.rot_bins_2[i].z();
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
            auto xr = loc.x() + radial_lines.rot_bins_3[i].x();
            auto yr = loc.y() + radial_lines.rot_bins_3[i].y();
            auto zr = loc.z() + radial_lines.rot_bins_3[i].z();
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

template<bool detect_surface, bool unity_width>
GridExcludedVolume helper(observer_ptr<grid::Grid> grid) {
    GridExcludedVolume vol;
    vol.interior.reserve(grid->get_volume());

    grid::detail::State exv_state;
    switch (settings::grid::exv::expansion_strategy) {
        case settings::grid::exv::ExvType::AtomicOnly:
            exv_state = grid::detail::State::VOLUME   | 
                        grid::detail::State::A_CENTER | 
                        grid::detail::State::A_AREA;
            break;

        case settings::grid::exv::ExvType::AtomicAndWater:
            exv_state = grid::detail::State::VOLUME   | 
                        grid::detail::State::A_CENTER | 
                        grid::detail::State::A_AREA   | 
                        grid::detail::State::W_CENTER | 
                        grid::detail::State::W_AREA;
            break;

        default:
            throw std::runtime_error("RawGridExv: Unknown expansion strategy. Did you forget to add a case?");
    }

    int stride = std::max(1., std::round(settings::grid::exv::width/settings::grid::cell_width));
    int buffer = std::max(1., std::round(std::max(settings::grid::min_exv_radius, 2.)/settings::grid::cell_width));

    const auto& axes = grid->get_axes();
    auto& gobj = grid->grid;
    auto[vmin, vmax] = grid->bounding_box_index();
    for (int i = std::max<int>(vmin.x()-buffer, 0); i < std::min<int>(vmax.x()+buffer+1, axes.x.bins); i+=stride) {
        for (int j = std::max<int>(vmin.y()-buffer, 0); j < std::min<int>(vmax.y()+buffer+1, axes.y.bins); j+=stride) {
            for (int k = std::max<int>(vmin.z()-buffer, 0); k < std::min<int>(vmax.z()+buffer+1, axes.z.bins); k+=stride) {
                auto& val = gobj.index(i, j, k);
                if (!(val & exv_state)) {continue;}

                // if we're not detecting the surface, everything is interior
                if constexpr (!detect_surface) {
                    vol.interior.emplace_back(grid->to_xyz(i, j, k));
                    continue;
                }

                // with non-unity widths we don't need the actual vectors yet since we have more work to do later
                auto collision = collision_check(grid, {i, j, k});
                if constexpr (!unity_width) {
                    if (!collision) {val |= grid::detail::RESERVED_1;}
                } else { // with unity widths our work is already done here
                    if (collision) {vol.interior.emplace_back(grid->to_xyz(i, j, k));}
                    else           {vol.surface.emplace_back( grid->to_xyz(i, j, k));}
                }
            }
        }
    }

    if constexpr (!unity_width) {
        int expand = std::round(settings::grid::exv::surface_thickness/settings::grid::cell_width)*0.5;
        int expand2 = expand*expand;
        auto mark_adjacent = [&gobj, exv_state, expand, expand2] (int i, int j, int k) {
            Vector3<int> origin{i, j, k};
            for (int l = -expand; l <= expand; l++) {
                for (int m = -expand; m <= expand; m++) {
                    for (int n = -expand; n <= expand; n++) {
                        auto& val = gobj.index(i+l, j+m, k+n);
                        if ((val & exv_state) && (origin.distance2(Vector3<int>{i+l, j+m, k+n}) <= expand2)) {
                            val |= grid::detail::RESERVED_2;
                        }
                    }
                }
            }
        };

        // expand the area around each surface voxel
        for (int i = std::max<int>(vmin.x()-buffer, 0); i < std::min<int>(vmax.x()+buffer+1, axes.x.bins); i+=stride) {
            for (int j = std::max<int>(vmin.y()-buffer, 0); j < std::min<int>(vmax.y()+buffer+1, axes.y.bins); j+=stride) {
                for (int k = std::max<int>(vmin.z()-buffer, 0); k < std::min<int>(vmax.z()+buffer+1, axes.z.bins); k+=stride) {
                    if (gobj.index(i, j, k) & grid::detail::RESERVED_1) {
                        mark_adjacent(i, j, k);
                    }
                }
            }
        }

        // collect the surface voxels
        for (int i = std::max<int>(vmin.x()-buffer, 0); i < std::min<int>(vmax.x()+buffer+1, axes.x.bins); i+=stride) {
            for (int j = std::max<int>(vmin.y()-buffer, 0); j < std::min<int>(vmax.y()+buffer+1, axes.y.bins); j+=stride) {
                for (int k = std::max<int>(vmin.z()-buffer, 0); k < std::min<int>(vmax.z()+buffer+1, axes.z.bins); k+=stride) {
                    auto& val = gobj.index(i, j, k);
                    if (val & (grid::detail::RESERVED_1 | grid::detail::RESERVED_2)) {
                        vol.surface.emplace_back(grid->to_xyz(i, j, k));
                        val &= ~(grid::detail::RESERVED_1 | grid::detail::RESERVED_2);
                        continue;
                    }

                    if (!(val & exv_state)) {continue;}
                    vol.interior.emplace_back(grid->to_xyz(i, j, k));
                }
            }
        }
    }

    assert(
        static_cast<int>(vol.interior.size() + vol.surface.size()) == grid->get_volume_bins() 
        && "RawGridExv: The number of interior and surface atoms does not match the number of volume bins."
    );
    logging::log(
        "RawGridWithSurfaceExv::create: added " + std::to_string(vol.interior.size()) + std::to_string(vol.surface.size()) + 
        "/" + std::to_string(grid->get_volume_bins()) + " atoms to the excluded volume."
        "\n\tOf these, " + std::to_string(vol.interior.size()) + " are interior atoms, and " + std::to_string(vol.surface.size()) + " are surface atoms."
    );

    return vol;
}

GridExcludedVolume RawGridWithSurfaceExv::create(observer_ptr<grid::Grid> grid, bool detect_surface) {
    int expand = std::round(settings::grid::exv::surface_thickness/settings::grid::cell_width);
    if (expand != 1) {
        return detect_surface ? helper<true, false>(grid) : helper<false, false>(grid);
    }
    return detect_surface ? helper<true, true>(grid) : helper<false, true>(grid);
}