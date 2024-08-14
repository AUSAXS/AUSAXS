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
#include <settings/MoleculeSettings.h>

#include <cassert>
#include <random>

using namespace data::record;

std::function<Vector3<double>()> hydrate::RadialHydration::noise_generator = [] () {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<> gauss(0, 0.75);
    return Vector3<double>(gauss(gen), gauss(gen), gauss(gen));
};

hydrate::RadialHydration::RadialHydration(observer_ptr<data::Molecule> protein) : GridBasedHydration(protein) {
    initialize();
    prepare_rotations();
}

hydrate::RadialHydration::RadialHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy) : GridBasedHydration(protein, std::move(culling_strategy)) {
    initialize();
    prepare_rotations();
}

hydrate::RadialHydration::~RadialHydration() = default;

void hydrate::RadialHydration::initialize() {
    hydrate::GridBasedHydration::initialize();
}

void hydrate::RadialHydration::set_noise_generator(std::function<Vector3<double>()>&& f) {
    noise_generator = f;
}

std::vector<grid::GridMember<data::record::Water>> hydrate::RadialHydration::generate_explicit_hydration() {
    assert(protein != nullptr && "RadialHydration::generate_explicit_hydration: protein is nullptr.");
    auto grid = protein->get_grid();
    assert(grid != nullptr && "RadialHydration::generate_explicit_hydration: grid is nullptr.");

    // we define a helper lambda
    std::vector<grid::GridMember<Water>> placed_water; 
    placed_water.reserve(grid->a_members.size());
    auto add_loc = [&] (Vector3<double>&& exact_loc) {
        Water a = Water::create_new_water(std::move(exact_loc));
        grid::GridMember<Water> gm = grid->add(a, true);
        placed_water.emplace_back(std::move(gm));
    };

    double rh = grid->get_hydration_radius() + settings::hydrate::shell_correction;
    for (const auto& atom : grid->a_members) {
        const auto& coords_abs = atom.get_atom().get_coordinates();
        double ra = grid->get_atomic_radius(atom.get_atom_type());
        double reff = ra + rh;
    
        for (unsigned int i = 0; i < rot_locs.size(); i++) {
            auto noise = noise_generator();
            auto bins = grid->to_bins_bounded(coords_abs + rot_locs[i]*reff + noise);
            if (grid->grid.is_empty_or_volume(bins.x(), bins.y(), bins.z()) && collision_check(Vector3<int>(bins.x(), bins.y(), bins.z()))) {
                Vector3<double> exact_loc = atom.get_atom().get_coordinates() + rot_locs[i]*reff + noise;
                add_loc(std::move(exact_loc));
            }
        }
    }
    return placed_water;
}

void hydrate::RadialHydration::prepare_rotations(int divisions) {
    auto grid = protein->get_grid();
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
            double x = std::cos(phi)*std::sin(theta);
            double y = std::sin(phi)*std::sin(theta);
            double z = std::cos(theta);
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

    double rh = grid->get_hydration_radius() + settings::hydrate::shell_correction;
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

bool hydrate::RadialHydration::collision_check(const Vector3<int>& loc) const {
    auto grid = protein->get_grid();
    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    int score = 0;

    // check if a location is out-of-bounds
    auto is_out_of_bounds = [&bins] (Vector3<int> v) {
        if (v.x() < 0 || bins.x() <= v.x() ) {return true;}
        if (v.y() < 0 || bins.y() <= v.y() ) {return true;}
        if (v.z() < 0 || bins.z() <= v.z() ) {return true;}
        return false;
    };

    unsigned int inside_1rh = 0;
    for (unsigned int i = 0; i < rot_locs.size(); i++) {
        {   // check for collisions at 1rh
            int xr = loc.x() + rot_bins_1rh[i].x();
            int yr = loc.y() + rot_bins_1rh[i].y();
            int zr = loc.z() + rot_bins_1rh[i].z();

            std::clamp(xr, 0, bins.x()-1);
            std::clamp(yr, 0, bins.y()-1);
            std::clamp(zr, 0, bins.z()-1);

            if (!gref.is_empty_or_volume(xr, yr, zr)) {
                if (2 < ++inside_1rh) {
                    return false;
                }
                continue;
            }
        }
        {   // check for collisions at 3rh
            int xr = loc.x() + rot_bins_3rh[i].x();
            int yr = loc.y() + rot_bins_3rh[i].y();
            int zr = loc.z() + rot_bins_3rh[i].z();
            if (is_out_of_bounds({xr, yr, zr})) {
                score += 2;
                continue;
            }

            if (!gref.is_empty_or_volume(xr, yr, zr)) {
                score -= 2;
                continue;
            }
            score += 2;
        }
        {   // check for collisions at 5rh
            int xr = loc.x() + rot_bins_5rh[i].x();
            int yr = loc.y() + rot_bins_5rh[i].y();
            int zr = loc.z() + rot_bins_5rh[i].z();
            if (is_out_of_bounds({xr, yr, zr})) {
                score += 1;
                continue;
            }

            if (!gref.is_empty_or_volume(xr, yr, zr)) {
                score -= 1;
                continue;
            }
        }
    }

    double max_points = 3*rot_bins_1rh.size();
    if (score <= settings::grid::detail::min_score*max_points) {
        return false;
    }
    return true;
}