#include <hydrate/placement/AxesPlacement.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <data/record/Water.h>
#include <math/Vector3.h>
#include <constants/Constants.h>

using namespace data::record;

#include <sstream>
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

    static int counter;
    counter = 0;
    // loop over the location of all member atoms
    double rh = grid->get_hydration_radius(); // radius of a water molecule
    std::cout << "(AxesPlacement) Grid[23, 32, 54] = " << (int) grid->grid.index(23, 32, 54) << std::endl;
    for (const auto& atom : grid->a_members) {
        double ra = grid->get_atomic_radius(atom.get_atom_type()); // radius of the atom
        double r_eff_real = ra+rh; // the effective bin radius
        int r_eff_bin = std::round(r_eff_real)/grid->get_width(); // the effective bin radius in bins

        auto coords_abs = atom.get_atom().get_coordinates();
        auto x = atom.get_bin_loc().x(), y = atom.get_bin_loc().y(), z = atom.get_bin_loc().z();
        auto xx = coords_abs.x(), yx = coords_abs.y(), zx = coords_abs.z(); // eXact coordinates
        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
        int xm = std::max<int>(std::round(xx-r_eff_bin), 0), xp = std::min<int>(std::round(xx+r_eff_bin), bins[0]-1); // xminus and xplus
        int ym = std::max<int>(std::round(yx-r_eff_bin), 0), yp = std::min<int>(std::round(yx+r_eff_bin), bins[1]-1); // yminus and yplus
        int zm = std::max<int>(std::round(zx-r_eff_bin), 0), zp = std::min<int>(std::round(zx+r_eff_bin), bins[2]-1); // zminus and zplus

        // check collisions for x ± r_eff
        std::stringstream ss;
        std::cout << "\nChecking x: " << std::endl;
        if ((gref.index(xm, y, z) == detail::EMPTY) && collision_check(Vector3<unsigned int>(xm, y, z), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.x() -= r_eff_real;
            std::cout << "\tAdding water at bin (" << xm << ", " << y << ", " << z << "), exact: " << exact_loc << std::endl;
            counter++;
            add_loc(exact_loc);
        } else if (gref.index(xm, y, z) != detail::EMPTY) {
            ss << "\tBin location (" << xm << ", " << y << ", " << z << ") already taken by " << (int) gref.index(xm, y, z) << std::endl;
        }
        if ((gref.index(xp, y, z) == detail::EMPTY) && collision_check(Vector3<unsigned int>(xp, y, z), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.x() += r_eff_real;
            std::cout << "\tAdding water at bin (" << xp << ", " << y << ", " << z << "), exact: " << exact_loc << std::endl;
            counter++;
            add_loc(exact_loc);
        } else if (gref.index(xp, y, z) != detail::EMPTY) {
            ss << "\tBin location (" << xp << ", " << y << ", " << z << ") already taken by " << (int) gref.index(xp, y, z) << std::endl;
        }

        // check collisions for y ± r_eff
        std::cout << "Checking y: " << std::endl;
        if ((gref.index(x, ym, z) == detail::EMPTY) && collision_check(Vector3<unsigned int>(x, ym, z), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.y() -= r_eff_real;
            std::cout << "\tAdding water at bin (" << x << ", " << ym << ", " << z << "), exact: " << exact_loc << std::endl;
            counter++;
            add_loc(exact_loc);
        } else if (gref.index(x, ym, z) != detail::EMPTY) {
            ss << "\tBin location (" << x << ", " << ym << ", " << z << ") already taken by " << (int) gref.index(x, ym, z) << std::endl;
        }

        if ((gref.index(x, yp, z) == detail::EMPTY) && collision_check(Vector3<unsigned int>(x, yp, z), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.y() += r_eff_real;
            std::cout << "\tAdding water at bin (" << x << ", " << yp << ", " << z << "), exact: " << exact_loc << std::endl;
            counter++;
            add_loc(exact_loc);
        } else if (gref.index(x, yp, z) != detail::EMPTY) {
            ss << "\tBin location (" << x << ", " << yp << ", " << z << ") already taken by " << (int) gref.index(x, yp, z) << std::endl;
        }

        // check collisions for z ± r_eff
        std::cout << "Checking z: " << std::endl;
        if ((gref.index(x, y, zm) == detail::EMPTY) && collision_check(Vector3<unsigned int>(x, y, zm), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.z() -= r_eff_real;
            std::cout << "\tAdding water at bin (" << x << ", " << y << ", " << zm << "), exact: " << exact_loc << std::endl;
            counter++;
            add_loc(exact_loc);
        } else if (gref.index(x, y, zm) != detail::EMPTY) {
            ss << "\tBin location (" << x << ", " << y << ", " << zm << ") already taken by " << (int) gref.index(x, y, zm) << std::endl;
        }

        if ((gref.index(x, y, zp) == detail::EMPTY) && collision_check(Vector3<unsigned int>(x, y, zp), ra)) {
            Vector3 exact_loc = atom.get_atom().get_coordinates();
            exact_loc.z() += r_eff_real;
            std::cout << "\tAdding water at bin (" << x << ", " << y << ", " << zp << "), exact: " << exact_loc << std::endl;
            counter++;
            add_loc(exact_loc);
        } else if (gref.index(x, y, zp) != detail::EMPTY) {
            ss << "\tBin location (" << x << ", " << y << ", " << zp << ") already taken by " << (int) gref.index(x, y, zp) << std::endl;
        }

        if (1560 <= counter && counter <= 1580) {
            std::cout << ss.str() << std::endl;
        }
    }
    std::cout << "(AxesPlacement) Grid[23, 32, 54] = " << (int) grid->grid.index(23, 32, 54) << std::endl;

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