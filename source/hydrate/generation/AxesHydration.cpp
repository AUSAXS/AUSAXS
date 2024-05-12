/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/generation/AxesHydration.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <math/Vector3.h>
#include <constants/Constants.h>

using namespace data::record;

hydrate::AxesHydration::AxesHydration(observer_ptr<data::Molecule> protein) : GridBasedHydration(protein) {
    initialize();
}

hydrate::AxesHydration::~AxesHydration() = default;

void hydrate::AxesHydration::initialize() {
    hydrate::GridBasedHydration::initialize();
    grid = protein->get_grid();
}

std::vector<grid::GridMember<data::record::Water>> hydrate::AxesHydration::generate_explicit_hydration() {
    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    // short lambda to actually place the generated water molecules
    std::vector<grid::GridMember<data::record::Water>> placed_water;
    placed_water.reserve(grid->a_members.size());
    auto add_loc = [&] (Vector3<double> exact_loc) {
        Water a = Water::create_new_water(exact_loc);
        grid::GridMember<Water> gm = grid->add(a, true);
        placed_water.emplace_back(std::move(gm));
    };

    // loop over the location of all member atoms
    double rh = grid->get_hydration_radius(); // radius of a water molecule
    for (const auto& atom : grid->a_members) {
        double ra = grid->get_atomic_radius(atom.get_atom_type()); // radius of the atom
        double r_eff_real = ra+rh; // the effective bin radius
        // int r_eff_bin = std::round(r_eff_real)/grid->get_width(); // the effective bin radius in bins

        auto coords_abs = atom.get_atom().get_coordinates();
        auto x = atom.get_bin_loc().x(), y = atom.get_bin_loc().y(), z = atom.get_bin_loc().z();

        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
        auto bin_min = grid->to_bins(coords_abs - r_eff_real);
        auto bin_max = grid->to_bins(coords_abs + r_eff_real);
        bin_min.x() = std::max<int>(bin_min.x(), 0); bin_max.x() = std::min<int>(bin_max.x(), bins[0]-1);
        bin_min.y() = std::max<int>(bin_min.y(), 0); bin_max.y() = std::min<int>(bin_max.y(), bins[1]-1);
        bin_min.z() = std::max<int>(bin_min.z(), 0); bin_max.z() = std::min<int>(bin_max.z(), bins[2]-1);

        // check collisions for x ± r_eff
        if ((gref.is_empty(bin_min.x(), y, z)) && collision_check(Vector3<unsigned int>(bin_min.x(), y, z), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.x() -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref.is_empty(bin_max.x(), y, z)) && collision_check(Vector3<unsigned int>(bin_max.x(), y, z), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.x() += r_eff_real;
            add_loc(exact_loc);
        }

        // check collisions for y ± r_eff
        if ((gref.is_empty(x, bin_min.y(), z)) && collision_check(Vector3<unsigned int>(x, bin_min.y(), z), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.y() -= r_eff_real;
            add_loc(exact_loc);
        }

        if ((gref.is_empty(x, bin_max.y(), z)) && collision_check(Vector3<unsigned int>(x, bin_max.x(), z), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.y() += r_eff_real;
            add_loc(exact_loc);
        }

        // check collisions for z ± r_eff
        if ((gref.is_empty(x, y, bin_min.z())) && collision_check(Vector3<unsigned int>(x, y, bin_min.z()), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.z() -= r_eff_real;
            add_loc(exact_loc);
        }

        if ((gref.is_empty(x, y, bin_max.z())) && collision_check(Vector3<unsigned int>(x, y, bin_max.z()), ra)) {
            Vector3 exact_loc = coords_abs;
            exact_loc.z() += r_eff_real;
            add_loc(exact_loc);
        }
    }
    return placed_water;
}

bool hydrate::AxesHydration::collision_check(const Vector3<unsigned int>& loc, double ra) const {
    static double rh = constants::radius::get_vdw_radius(constants::atom_t::O); // radius of a water molecule
    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    
    int x = loc.x(), y = loc.y(), z = loc.z();

    // loop over the box [x-r, x+r][y-r, y+r][z-r, z+r]
    int r = gref.is_atom_center(x, y, z) ? ra : rh;

    // we use the range (x-r) to (x+r+1) since the first is inclusive and the second is exclusive. 
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) bins[0])-1; // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) bins[1])-1; // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) bins[2])-1; // zminus and zplus
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                if (!gref.is_empty(i, j, k) && pow(x-i, 2) + pow(y-j, 2) + pow(z-k, 2) < r*r) {return false;}
            }
        }
    }
    return true;
}