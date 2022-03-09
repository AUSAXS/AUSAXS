#include "hydrate/AxesPlacement.h"
#include "hydrate/Grid.h"

vector<grid::GridMember<Hetatom>> grid::AxesPlacement::place() const {
    // dereference the values we'll need for better performance
    vector<vector<vector<char>>>& gref = grid->grid;
    const vector<int> bins = grid->get_bins();
    const int ra = grid->ra; const int rh = grid->rh;

    // short lambda to actually place the generated water molecules
    vector<GridMember<Hetatom>> placed_water(grid->a_members.size());
    size_t index = 0;
    auto add_loc = [&] (const Vector3 exact_loc) {
        Hetatom a = Hetatom::create_new_water(exact_loc);
        GridMember<Hetatom> gm = grid->add(a, true);
        if (__builtin_expect(placed_water.size() <= index, false)) {
            placed_water.resize(2*index);
        }
        placed_water[index++] = gm;
    };

    // loop over the location of all member atoms
    int r_eff = ra+rh;                           // the effective bin radius
    double r_eff_real = r_eff*grid->get_width(); // the effective real radius
    for (const auto& atom : grid->a_members) {
        const int x = atom.loc[0], y = atom.loc[1], z = atom.loc[2];

        // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
        int xm = std::max(x-r_eff, 0), xp = std::min(x+r_eff, bins[0]-1); // xminus and xplus
        int ym = std::max(y-r_eff, 0), yp = std::min(y+r_eff, bins[1]-1); // yminus and yplus
        int zm = std::max(z-r_eff, 0), zp = std::min(z+r_eff, bins[2]-1); // zminus and zplus

        // check collisions for x ± r_eff
        if ((gref[xm][y][z] == 0) && collision_check({xm, y, z})) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.x -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref[xp][y][z] == 0) && collision_check({xp, y, z})) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.x += r_eff_real;
            add_loc(exact_loc);
        }

        // check collisions for y ± r_eff
        if ((gref[x][ym][z] == 0) && collision_check({x, ym, z})) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.y -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref[x][yp][z] == 0) && collision_check({x, yp, z})) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.y += r_eff_real;
            add_loc(exact_loc);
        }

        // check collisions for z ± r_eff
        if ((gref[x][y][zm] == 0) && collision_check({x, y, zm})) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.z -= r_eff_real;
            add_loc(exact_loc);
        }
        if ((gref[x][y][zp] == 0) && collision_check({x, y, zp})) {
            Vector3 exact_loc = atom.atom.coords;
            exact_loc.z += r_eff_real;
            add_loc(exact_loc);
        }
    }

    placed_water.resize(index);
    return placed_water;
}

bool grid::AxesPlacement::collision_check(const vector<int> loc) const {
    // dereference the values we'll need for better performance
    vector<vector<vector<char>>>& gref = grid->grid;
    const vector<int> bins = grid->get_bins();
    const int ra = grid->ra, rh = grid->rh;
    const int x = loc[0], y = loc[1], z = loc[2];

    // loop over the box [x-r, x+r][y-r, y+r][z-r, z+r]
    int r = gref[x][y][z] == 'A' ? ra : rh;

    // we use the range (x-r) to (x+r+1) since the first is inclusive and the second is exclusive. 
    int xm = std::max(x-r, 0), xp = std::min(x+r+1, bins[0]-1); // xminus and xplus
    int ym = std::max(y-r, 0), yp = std::min(y+r+1, bins[1]-1); // yminus and yplus
    int zm = std::max(z-r, 0), zp = std::min(z+r+1, bins[2]-1); // zminus and zplus
    for (int i = xm; i < xp; i++) {
        for (int j = ym; j < yp; j++) {
            for (int k = zm; k < zp; k++) {
                if (gref[i][j][k] != 0 && sqrt(pow(x-i, 2) + pow(y-j, 2) + pow(z-k, 2)) < r) {return false;}
            }
        }
    }
    return true;
}