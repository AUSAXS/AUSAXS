#include <hydrate/placement/AxesPlacement.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <data/record/Water.h>
#include <math/Vector3.h>
#include <constants/Constants.h>

using namespace data::record;

std::vector<grid::GridMember<data::record::Water>> grid::AxesPlacement::place() const {
    detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    // short lambda to actually place the generated water molecules
    std::vector<GridMember<data::record::Water>> placed_water(grid->a_members.size());
    unsigned int index = 0;
    auto add_loc = [&] (Vector3<double> exact_loc) {
        Water a = Water::create_new_water(exact_loc);
        GridMember<Water> gm = grid->add(a, true);
        if (placed_water.size() <= index) [[unlikely]] {
            placed_water.resize(2*index);
        }
        placed_water[index++] = gm;
    };

    // loop over the location of all member atoms
    double rh = grid->get_hydration_radius(); // radius of a water molecule
    for (const auto& atom : grid->a_members) {
        int x = atom.get_bin_loc().x(), y = atom.get_bin_loc().y(), z = atom.get_bin_loc().z();
        double ra = grid->get_atomic_radius(atom.get_atom_type()); // radius of the atom
        double r_eff_real = ra+rh; // the effective bin radius
        int r_eff_bin = std::round(r_eff_real)/grid->get_width(); // the effective bin radius in bins

        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
        unsigned int xm = std::max(x-r_eff_bin, 0), xp = std::min(x+r_eff_bin, (int) bins[0]-1); // xminus and xplus
        unsigned int ym = std::max(y-r_eff_bin, 0), yp = std::min(y+r_eff_bin, (int) bins[1]-1); // yminus and yplus
        unsigned int zm = std::max(z-r_eff_bin, 0), zp = std::min(z+r_eff_bin, (int) bins[2]-1); // zminus and zplus

        // check collisions for x ± r_eff
        if ((gref.index(xm, y, z) == detail::EMPTY) && collision_check(Vector3<unsigned int>(xm, y, z), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.x() -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref.index(xp, y, z) == detail::EMPTY) && collision_check(Vector3<unsigned int>(xp, y, z), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.x() += r_eff_real;
            add_loc(exact_loc);
        }

        // check collisions for y ± r_eff
        if ((gref.index(x, ym, z) == detail::EMPTY) && collision_check(Vector3<unsigned int>(x, ym, z), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.y() -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref.index(x, yp, z) == detail::EMPTY) && collision_check(Vector3<unsigned int>(x, yp, z), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.y() += r_eff_real;
            add_loc(exact_loc);
        }

        // check collisions for z ± r_eff
        if ((gref.index(x, y, zm) == detail::EMPTY) && collision_check(Vector3<unsigned int>(x, y, zm), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.z() -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref.index(x, y, zp) == detail::EMPTY) && collision_check(Vector3<unsigned int>(x, y, zp), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.z() += r_eff_real;
            add_loc(exact_loc);
        }
    }

    placed_water.resize(index);
    return placed_water;
}

bool grid::AxesPlacement::collision_check(const Vector3<unsigned int>& loc, double ra) const {
    static double rh = constants::radius::get_vdw_radius(constants::atom_t::O); // radius of a water molecule

    detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    
    int x = loc.x(), y = loc.y(), z = loc.z();

    // loop over the box [x-r, x+r][y-r, y+r][z-r, z+r]
    int r = gref.index(x, y, z) == detail::A_CENTER ? ra : rh;

    // we use the range (x-r) to (x+r+1) since the first is inclusive and the second is exclusive. 
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) bins[0])-1; // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) bins[1])-1; // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) bins[2])-1; // zminus and zplus
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                if (gref.index(i, j, k) != detail::EMPTY && pow(x-i, 2) + pow(y-j, 2) + pow(z-k, 2) < r*r) {return false;}
            }
        }
    }
    return true;
}