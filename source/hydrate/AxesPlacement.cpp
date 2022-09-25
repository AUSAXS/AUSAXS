#include <hydrate/AxesPlacement.h>
#include <hydrate/Grid.h>

using std::vector;

vector<grid::GridMember<Hetatom>> grid::AxesPlacement::place() const {
    // dereference the values we'll need for better performance
    GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    unsigned int ra = grid->ra, rh = grid->rh;

    // short lambda to actually place the generated water molecules
    vector<GridMember<Hetatom>> placed_water(grid->a_members.size());
    size_t index = 0;
    auto add_loc = [&] (Vector3<double> exact_loc) {
        Hetatom a = Hetatom::create_new_water(exact_loc);
        GridMember<Hetatom> gm = grid->add(a, true);
        if (__builtin_expect(placed_water.size() <= index, false)) {
            placed_water.resize(2*index);
        }
        placed_water[index++] = gm;
    };

    // loop over the location of all member atoms
    int r_eff = ra+rh;                  // the effective bin radius
    double r_eff_real = r_eff*grid->get_width(); // the effective real radius
    for (const auto& atom : grid->a_members) {
        int x = atom.loc.x(), y = atom.loc.y(), z = atom.loc.z();

        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
        unsigned int xm = std::max(x-r_eff, 0), xp = std::min(x+r_eff, (int) bins[0]-1); // xminus and xplus
        unsigned int ym = std::max(y-r_eff, 0), yp = std::min(y+r_eff, (int) bins[1]-1); // yminus and yplus
        unsigned int zm = std::max(z-r_eff, 0), zp = std::min(z+r_eff, (int) bins[2]-1); // zminus and zplus

        // check collisions for x ± r_eff
        if ((gref.index(xm, y, z) == GridObj::EMPTY) && collision_check(Vector3<unsigned int>(xm, y, z))) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.x() -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref.index(xp, y, z) == GridObj::EMPTY) && collision_check(Vector3<unsigned int>(xp, y, z))) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.x() += r_eff_real;
            add_loc(exact_loc);
        }

        // check collisions for y ± r_eff
        if ((gref.index(x, ym, z) == GridObj::EMPTY) && collision_check(Vector3<unsigned int>(x, ym, z))) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.y() -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref.index(x, yp, z) == GridObj::EMPTY) && collision_check(Vector3<unsigned int>(x, yp, z))) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.y() += r_eff_real;
            add_loc(exact_loc);
        }

        // check collisions for z ± r_eff
        if ((gref.index(x, y, zm) == GridObj::EMPTY) && collision_check(Vector3<unsigned int>(x, y, zm))) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.z() -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref.index(x, y, zp) == GridObj::EMPTY) && collision_check(Vector3<unsigned int>(x, y, zp))) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.z() += r_eff_real;
            add_loc(exact_loc);
        }
    }

    placed_water.resize(index);
    return placed_water;
}

bool grid::AxesPlacement::collision_check(const Vector3<unsigned int>& loc) const {
    // dereference the values we'll need for better performance
    GridObj& gref = grid->grid;
    auto bins = grid->get_bins();
    int ra = grid->ra, rh = grid->rh;
    int x = loc.x(), y = loc.y(), z = loc.z();

    // loop over the box [x-r, x+r][y-r, y+r][z-r, z+r]
    int r = gref.index(x, y, z) == GridObj::A_CENTER ? ra : rh;

    // we use the range (x-r) to (x+r+1) since the first is inclusive and the second is exclusive. 
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, (int) bins[0])-1; // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, (int) bins[1])-1; // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, (int) bins[2])-1; // zminus and zplus
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                if (gref.index(i, j, k) != GridObj::EMPTY && sqrt(pow(x-i, 2) + pow(y-j, 2) + pow(z-k, 2)) < r) {return false;}
            }
        }
    }
    return true;
}