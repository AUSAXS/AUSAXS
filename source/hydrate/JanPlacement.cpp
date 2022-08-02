#include <hydrate/JanPlacement.h>
#include <hydrate/Grid.h>

using std::vector;

vector<grid::GridMember<Hetatom>> grid::JanPlacement::place() const {
    // dereference the values we'll need for better performance
    vector<vector<vector<char>>>& gref = grid->grid;
    auto bins = grid->get_bins();

    // place a water molecule (note: not added to the grid before the end of this method)
    vector<Hetatom> placed_water(grid->a_members.size());
    size_t index = 0;
    auto add_loc = [&] (const vector<unsigned int> v) {
        Hetatom a = Hetatom::create_new_water(grid->to_xyz(v));
        if (__builtin_expect(placed_water.size() <= index, false)) {
            placed_water.resize(2*index);
        }
        placed_water[index++] = a;
    };

    // loop over the location of all member atoms
    int r_eff = grid->ra;
    auto[min, max] = grid->bounding_box_index();
    for (unsigned int i = min[0]; i < max[0]; i++) {
        for (unsigned int j = min[1]; j < max[1]; j++) {
            for (unsigned int k = min[2]; k < max[2]; k++) {
                if (gref[i][j][k] == 0) {continue;}

                // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
                unsigned int im = std::max(i-r_eff, 0U), ip = std::min(i+r_eff, bins[0]-1); // xminus and xplus
                unsigned int jm = std::max(j-r_eff, 0U), jp = std::min(j+r_eff, bins[1]-1); // yminus and yplus
                unsigned int km = std::max(k-r_eff, 0U), kp = std::min(k+r_eff, bins[2]-1); // zminus and zplus

                // check collisions for x ± r_eff
                if (gref[im][j][k] == 0) {add_loc({im, j, k});}
                if (gref[ip][j][k] == 0) {add_loc({ip, j, k});}

                // check collisions for y ± r_eff
                if (gref[i][jm][k] == 0) {add_loc({i, jm, k});}
                if (gref[i][jp][k] == 0) {add_loc({i, jp, k});}

                // check collisions for z ± r_eff
                if (gref[i][j][km] == 0) {add_loc({i, j, km});}
                if (gref[i][j][kp] == 0) {add_loc({i, j, kp});}
            }
        }
    }

    // finally we can add the atoms to the grid
    placed_water.resize(index);
    vector<grid::GridMember<Hetatom>> v = grid->add(placed_water);
    grid->expand_volume();

    return v;
}