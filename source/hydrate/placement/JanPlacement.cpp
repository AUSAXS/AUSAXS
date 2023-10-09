#include <hydrate/placement/JanPlacement.h>
#include <hydrate/GridMember.h>
#include <hydrate/Grid.h>
#include <data/Water.h>
#include <utility/Constants.h>

std::vector<grid::GridMember<Water>> grid::JanPlacement::place() const {
    // dereference the values we'll need for better performance
    GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    // place a water molecule (note: not added to the grid before the end of this method)
    std::vector<Water> placed_water(grid->a_members.size());
    unsigned int index = 0;
    auto add_loc = [&] (const Vector3<int>& v) {
        Water a = Water::create_new_water(grid->to_xyz(v));
        if (placed_water.size() <= index) [[unlikely]] {
            placed_water.resize(2*index);
        }
        placed_water[index++] = a;
    };

    // loop over the location of all member atoms
    int r_eff = (grid->get_atomic_radius(constants::atom_t::C) + grid->get_hydration_radius())/grid->get_width();
    auto[min, max] = grid->bounding_box_index();
    for (int i = min.x(); i < max.x(); i++) {
        for (int j = min.y(); j < max.y(); j++) {
            for (int k = min.z(); k < max.z(); k++) {
                if (gref.index(i, j, k) == GridObj::EMPTY) {continue;}

                // we define a small box of size [i-rh, i+rh][j-rh, j+rh][z-rh, z+rh]
                int im = std::max(i-r_eff, 0), ip = std::min(i+r_eff, (int) bins.x()-1); // xminus and xplus
                int jm = std::max(j-r_eff, 0), jp = std::min(j+r_eff, (int) bins.y()-1); // yminus and yplus
                int km = std::max(k-r_eff, 0), kp = std::min(k+r_eff, (int) bins.z()-1); // zminus and zplus

                // check collisions for x ± r_eff                
                if (gref.index(im, j, k) == GridObj::EMPTY) {add_loc(Vector3<int>(im, j, k));}
                if (gref.index(ip, j, k) == GridObj::EMPTY) {add_loc(Vector3<int>(ip, j, k));}

                // check collisions for y ± r_eff
                if (gref.index(i, jp, k) == GridObj::EMPTY) {add_loc(Vector3<int>(i, jp, k));}
                if (gref.index(i, jm, k) == GridObj::EMPTY) {add_loc(Vector3<int>(i, jm, k));}

                // check collisions for z ± r_eff
                if (gref.index(i, j, km) == GridObj::EMPTY) {add_loc(Vector3<int>(i, j, km));}
                if (gref.index(i, j, kp) == GridObj::EMPTY) {add_loc(Vector3<int>(i, j, kp));}
            }
        }
    }

    // finally we can add the atoms to the grid
    placed_water.resize(index);
    std::vector<grid::GridMember<Water>> v = grid->add(placed_water);
    grid->expand_volume();

    return v;
}