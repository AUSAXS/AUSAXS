/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/detail/GridSurfaceDetection.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/record/Atom.h>
#include <settings/GridSettings.h>

#include <unordered_map>

using namespace grid::detail;

GridSurfaceDetection::GridSurfaceDetection(observer_ptr<grid::Grid> grid) : RadialLineGenerator(grid), grid(grid) {}

GridSurfaceDetection::~GridSurfaceDetection() = default;

bool GridSurfaceDetection::collision_check(const Vector3<int>& loc) {
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

    for (unsigned int i = 0; i < rot_locs.size(); i++) {
        {   // check for collisions at 1r
            auto xr = loc.x() + rot_bins_1rh[i].x();
            auto yr = loc.y() + rot_bins_1rh[i].y();
            auto zr = loc.z() + rot_bins_1rh[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 1;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                continue;
            }
        }

        {   // check for collisions at 3r
            auto xr = loc.x() + rot_bins_3rh[i].x();
            auto yr = loc.y() + rot_bins_3rh[i].y();
            auto zr = loc.z() + rot_bins_3rh[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 1;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                continue;
            }
        }

        {   // check for collisions at 5r
            auto xr = loc.x() + rot_bins_5rh[i].x();
            auto yr = loc.y() + rot_bins_5rh[i].y();
            auto zr = loc.z() + rot_bins_5rh[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 1;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                continue;
            }
        }

        {   // check for collisions at 7r
            auto xr = loc.x() + rot_bins_7rh[i].x();
            auto yr = loc.y() + rot_bins_7rh[i].y();
            auto zr = loc.z() + rot_bins_7rh[i].z();
            if (is_out_of_bounds({xr, yr, zr})) [[unlikely]] {
                score += 1;
                continue;
            }

            if (!gref.is_empty_or_water(xr, yr, zr)) {
                continue;
            }
        }

        score += 1;
    }
    std::cout << "Score: " << score << std::endl;
    return score < 16;
}

GridExcludedVolume GridSurfaceDetection::detect_atoms() {
    GridExcludedVolume vol;
    vol.interior.reserve(grid->a_members.size());
    vol.surface.reserve(std::pow(grid->a_members.size(), 2./3));

    for (const auto& atom : grid->a_members) {
        auto coords_abs = atom.get_atom().get_coordinates();
        auto coords_bin = grid->to_bins_bounded(coords_abs);
        if (collision_check(coords_bin)) {
            vol.surface.emplace_back(coords_abs);
        } else {
            vol.interior.emplace_back(coords_abs);
        }
    }

    return vol;
}

GridExcludedVolume GridSurfaceDetection::detect_voxels() {
    GridExcludedVolume vol;
    vol.interior.reserve(grid->get_volume());
    vol.surface.reserve(std::pow(grid->get_volume(), 2./3));

    int buffer = std::round(2./settings::grid::width);
    int stride = std::round(2*settings::grid::exv_radius/settings::grid::width);
    const auto& axes = grid->get_axes();
    const auto& gobj = grid->grid;
    auto[vmin, vmax] = grid->bounding_box_index();
    for (int i = std::max<int>(vmin.x()-buffer, 0); i < std::min<int>(vmax.x()+buffer, axes.x.bins); i+=stride) {
        for (int j = std::max<int>(vmin.y()-buffer, 0); j < std::min<int>(vmax.y()+buffer, axes.y.bins); j+=stride) {
            for (int k = std::max<int>(vmin.z()-buffer, 0); k < std::min<int>(vmax.z()+buffer, axes.z.bins); k+=stride) {
                switch (gobj.index(i, j, k)) {
                    case detail::VOLUME:
                    case detail::A_AREA:
                    case detail::A_CENTER:
                        if (collision_check({i, j, k})) {
                            vol.surface.emplace_back(grid->to_xyz(i, j, k));
                        }
                    default:
                        break;
                }
            }
        }
    }

    return vol;
}

// // https://stackoverflow.com/a/42701911
// struct ArrayHasher {
//     std::size_t operator()(const std::array<int, 3>& a) const {
//         std::size_t h = 0;

//         for (auto e : a) {
//             h ^= std::hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2); 
//         }
//         return h;
//     }   
// };

// GridExcludedVolume GridSurfaceDetection::detect_voxels() {
//     GridExcludedVolume vol;
//     vol.interior.reserve(grid->get_volume());
//     vol.surface.reserve(std::pow(grid->get_volume(), 2./3));
//     auto atoms = detect_atoms();
//     auto cm = std::accumulate(atoms.interior.begin(), atoms.interior.end(), Vector3<double>());
//     cm = std::accumulate(atoms.surface.begin(), atoms.surface.end(), cm)/(atoms.interior.size() + atoms.surface.size());

//     auto cm_bins = grid->to_bins(cm);
//     auto oneA = 1./grid->get_width(); // one Angstrom
//     int buffer = std::round(2./settings::grid::width);
//     int stride = std::round(2*settings::grid::exv_radius/settings::grid::width);
//     const auto& axes = grid->get_axes();
//     auto& gobj = grid->grid;

//     std::unordered_map<std::array<int, 3>, detail::State, ArrayHasher> previous_state;
//     for (const auto& atom : atoms.surface) {
//         auto atom_bins = grid->to_bins(atom);
//         auto distance = atom_bins.distance2(cm_bins);

//         for (int i = std::max<int>(atom_bins.x()-buffer, 0); i < std::min<int>(atom_bins.x()+buffer, axes.x.bins); i+=stride) {
//             for (int j = std::max<int>(atom_bins.y()-buffer, 0); j < std::min<int>(atom_bins.y()+buffer, axes.y.bins); j+=stride) {
//                 for (int k = std::max<int>(atom_bins.z()-buffer, 0); k < std::min<int>(atom_bins.z()+buffer, axes.z.bins); k+=stride) {
//                     auto& val = gobj.index(i, j, k);
//                     switch (val) {
//                         case detail::VOLUME:
//                         case detail::A_AREA:
//                         case detail::A_CENTER:
//                             // if (distance < cm_bins.distance2(Vector3{i, j, k})) {
//                                 previous_state.emplace(std::array{i, j, k}, val);
//                                 val = detail::RESERVED;
//                             // }
//                         default:
//                             break;
//                     }
//                 }
//             }
//         }
//     }

//     auto[vmin, vmax] = grid->bounding_box_index();
//     for (int i = std::max<int>(vmin.x()-buffer, 0); i < std::min<int>(vmax.x()+buffer, axes.x.bins); i+=stride) {
//         for (int j = std::max<int>(vmin.y()-buffer, 0); j < std::min<int>(vmax.y()+buffer, axes.y.bins); j+=stride) {
//             for (int k = std::max<int>(vmin.z()-buffer, 0); k < std::min<int>(vmax.z()+buffer, axes.z.bins); k+=stride) {
//                 auto& val = gobj.index(i, j, k);
//                 switch (val) {
//                     case detail::VOLUME:
//                     case detail::A_AREA:
//                     case detail::A_CENTER: 
//                         vol.interior.emplace_back(grid->to_xyz(i, j, k));
//                         break;
//                     case detail::RESERVED:
//                         vol.surface.emplace_back(grid->to_xyz(i, j, k));
//                         val = previous_state.at({i, j, k});
//                     default:
//                         break;
//                 }
//             }
//         }
//     }
//     return vol;
// }