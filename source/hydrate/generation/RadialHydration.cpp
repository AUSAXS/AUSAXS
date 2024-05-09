/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/generation/RadialHydration.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <constants/Constants.h>
#include <settings/GridSettings.h>

using namespace data::record;

hydrate::RadialHydration::RadialHydration(observer_ptr<data::Molecule> protein) : GridBasedHydration(protein) {
    initialize();
    prepare_rotations();
}

hydrate::RadialHydration::~RadialHydration() = default;

void hydrate::RadialHydration::initialize() {
    hydrate::GridBasedHydration::initialize();
    grid = protein->get_grid();
}

void hydrate::RadialHydration::prepare_rotations(int divisions) {
    double width = grid->get_width();

    std::vector<Vector3<int>> bins_1rh;
    std::vector<Vector3<int>> bins_3rh;
    std::vector<Vector3<int>> bins_5rh;
    std::vector<Vector3<int>> bins_7rh;
    std::vector<Vector3<double>> locs;
    double ang = 2*constants::pi/divisions;

    // we generate one octant of a sphere, and then reflect it to generate the rest
    // we do this to ensure the sphere is symmetric. If we simply generate it all at once, floating-point errors moves some of the bins around
    std::vector<Vector3<double>> sphere;
    for (double theta = 0; theta <= constants::pi*0.5; theta+=ang) {
        for (double phi = 0; phi <= constants::pi*0.5; phi+=ang) {
            double x = cos(phi)*sin(theta);
            double y = sin(phi)*sin(theta);
            double z = cos(theta);
            sphere.push_back({x, y, z});
            sphere.push_back({-x, y, z});
            sphere.push_back({x, -y, z});
            sphere.push_back({-x, -y, z});
            sphere.push_back({x, y, -z});
            sphere.push_back({-x, y, -z});
            sphere.push_back({x, -y, -z});
            sphere.push_back({-x, -y, -z});
        }
    }

    // remove duplicates
    std::vector<Vector3<double>> rots;
    for (auto& p : sphere) {
        bool present = false;
        for (int i = 0; i < 3; i++) { // fix the easy floating point errors
            if (abs(p[i]) < 1e-5) {p[i] = 0;}
        }
        for (const auto& r : rots) { // go through all rotations and try to find a duplicate entry
            if (r.distance(p) < 1e-5) {
                present = true;
                break;
            }
        }
        if (!present) { // if the element was not already present
            rots.push_back(p); // add it
        }
    }

    double rh = constants::radius::get_vdw_radius(constants::atom_t::O);
    for (const auto& rot : rots) {
        double xr = rot.x(), yr = rot.y(), zr = rot.z();
        bins_1rh.push_back(Vector3<int>(std::round(  rh*xr)/width, std::round(  rh*yr)/width, std::round(  rh*zr)/width));
        bins_3rh.push_back(Vector3<int>(std::round(3*rh*xr)/width, std::round(3*rh*yr)/width, std::round(3*rh*zr)/width));
        bins_5rh.push_back(Vector3<int>(std::round(5*rh*xr)/width, std::round(5*rh*yr)/width, std::round(5*rh*zr)/width));
        bins_7rh.push_back(Vector3<int>(std::round(7*rh*xr)/width, std::round(7*rh*yr)/width, std::round(7*rh*zr)/width));
        locs.push_back(Vector3<double>(xr, yr, zr));
    }

    // set the member vectors
    rot_bins_1rh = std::move(bins_1rh);
    rot_bins_3rh = std::move(bins_3rh);
    rot_bins_5rh = std::move(bins_5rh);
    rot_bins_7rh = std::move(bins_7rh);
    rot_locs = std::move(locs);
}

std::vector<grid::GridMember<data::record::Water>> hydrate::RadialHydration::generate_explicit_hydration() {
    // we define a helper lambda
    std::vector<grid::GridMember<Water>> placed_water; 
    placed_water.reserve(grid->a_members.size());
    auto add_loc = [&] (Vector3<double> exact_loc) {
        Water a = Water::create_new_water(exact_loc);
        grid::GridMember<Water> gm = grid->add(a, true);
        placed_water.emplace_back(std::move(gm));
    };

    double rh = grid->get_hydration_radius();
    for (const auto& atom : grid->a_members) {
        auto coords_abs = atom.get_atom().get_coordinates();
        double ra = grid->get_atomic_radius(atom.get_atom_type());
        double reff = ra + rh;
        for (unsigned int i = 0; i < rot_locs.size(); i++) {
            auto bins = grid->to_bins_bounded(coords_abs + rot_locs[i]*reff);

            // we have to make sure we don't check the direction of the atom we are trying to place this water on
            Vector3<int> skip_bin(bins.x()-rot_bins_1rh[i].x(), bins.y()-rot_bins_1rh[i].y(), bins.z()-rot_bins_1rh[i].z());
            if (grid->grid.is_empty_or_volume(bins.x(), bins.y(), bins.z()) && collision_check(Vector3<int>(bins.x(), bins.y(), bins.z()), skip_bin)) {
                Vector3<double> exact_loc = atom.get_atom().get_coordinates() + rot_locs[i]*reff;
                add_loc(exact_loc);
            }
        }
    }
    return placed_water;
}

bool hydrate::RadialHydration::collision_check(const Vector3<int>& loc, const Vector3<int>& skip_bin) const {
    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    int score = 0;

    // check if a location is out-of-bounds
    auto is_out_of_bounds = [&bins] (Vector3<int> v) {
        if (v.x() < 0 || (int) bins.x() <= v.x() ) {return true;}
        if (v.y() < 0 || (int) bins.y() <= v.y() ) {return true;}
        if (v.z() < 0 || (int) bins.z() <= v.z() ) {return true;}
        return false;
    };

    for (unsigned int i = 0; i < rot_locs.size(); i++) {
        // check for collisions at 1rh
        int xr = loc.x() + rot_bins_1rh[i].x();
        int yr = loc.y() + rot_bins_1rh[i].y();
        int zr = loc.z() + rot_bins_1rh[i].z();

        // check for bounds
        if (xr < 0) xr = 0;
        if (xr >= (int) bins.x()) xr = bins.x()-1;
        if (yr < 0) yr = 0;
        if (yr >= (int) bins.y()) yr = bins.y()-1;
        if (zr < 0) zr = 0;
        if (zr >= (int) bins.z()) zr = bins.z()-1;

        if (!gref.is_empty_or_volume(xr, yr, zr)) {
            if (Vector3(xr, yr, zr) == skip_bin) {continue;} // skip the bin containing the atom we're trying to place this water molecule on

            return false;
        }

        // check if we're in a cavity
        for (unsigned int j = 0; j < rot_bins_3rh.size(); j++) {
            // check at 3r
            int xr = loc.x() + rot_bins_3rh[j].x();
            int yr = loc.y() + rot_bins_3rh[j].y();
            int zr = loc.z() + rot_bins_3rh[j].z();
            if (is_out_of_bounds({xr, yr, zr})) {   // if the line goes out of bounds, we know for sure it won't intersect anything
                score += 3;                         // so we add three points and move on to the next
                continue;
            }
            if (!gref.is_empty_or_volume(xr, yr, zr)) { // if the line intersects something at 3r, we don't check the other two points of the same line
                score -= 3;                             // but immediately subtract 3 points and move on to the next
                continue;
            } else {
                score++;
            } 

            // check at 5r
            xr = loc.x() + rot_bins_5rh[j].x();
            yr = loc.y() + rot_bins_5rh[j].y();
            zr = loc.z() + rot_bins_5rh[j].z();
            if (is_out_of_bounds({xr, yr, zr})) {
                score += 2;
                continue;
            }
            if (!gref.is_empty_or_volume(xr, yr, zr)) {
                score -= 2;
                continue;
            } else {
                score++;
            } 

            // check at 7r
            xr = loc.x() + rot_bins_7rh[j].x();
            yr = loc.y() + rot_bins_7rh[j].y();
            zr = loc.z() + rot_bins_7rh[j].z();
            if (is_out_of_bounds({xr, yr, zr})) {
                score += 1;
                continue;
            }
            if (!gref.is_empty_or_volume(xr, yr, zr)) {
                score -= 1;
                continue;
            } else {
                score++;
            }
        }
    }

    if (score <= settings::grid::detail::min_score*rot_bins_1rh.size()) {
        return false;
    }
    return true;
}