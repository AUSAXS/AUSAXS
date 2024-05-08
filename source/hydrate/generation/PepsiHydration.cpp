/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/generation/PepsiHydration.h>
#include <data/record/Water.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <settings/GridSettings.h>
#include <settings/MoleculeSettings.h>

using namespace hydrate;
using namespace data::record;

PepsiHydration::PepsiHydration(observer_ptr<data::Molecule> protein) : GridBasedHydration(protein) {
    initialize();
}

void PepsiHydration::initialize() {
    std::cout << "Initializing PepsiHydration" << std::endl;
    settings::hydrate::culling_strategy = settings::hydrate::CullingStrategy::NoStrategy;
    GridBasedHydration::initialize();
    grid = protein->get_grid();
}

PepsiHydration::~PepsiHydration() = default;

void PepsiHydration::modified_expand_volume(grid::GridMember<data::record::Atom>& atom) {
    if (atom.is_expanded()) {return;} // check if this location has already been expanded
    atom.set_expanded(true); // mark this location as expanded

    grid::detail::GridObj& gref = grid->grid;
    auto axes = grid->get_axes();

    double r = 3/settings::grid::width; // fixed radius of 3Å

    // create a box of size [x-r, x+r][y-r, y+r][z-r, z+r] within the bounds
    int x = atom.get_bin_loc().x(), y = atom.get_bin_loc().y(), z = atom.get_bin_loc().z(); 
    double rvdw = r/settings::grid::width;
    double rvdw2 = std::pow(rvdw, 2);

    int xm = std::max<int>(x - std::ceil(r), 0), xp = std::min<int>(x + std::ceil(r) + 1, axes.x.bins); // xminus and xplus
    int ym = std::max<int>(y - std::ceil(r), 0), yp = std::min<int>(y + std::ceil(r) + 1, axes.y.bins); // yminus and yplus
    int zm = std::max<int>(z - std::ceil(r), 0), zp = std::min<int>(z + std::ceil(r) + 1, axes.z.bins); // zminus and zplus

    // loop over each bin in the box
    int added_volume = 0;

    // i, j, k *must* be ints to avoid unsigned underflow
    for (int i = xm; i < xp; ++i) {
        double x2 = std::pow(x - i, 2);
        for (int j = ym; j < yp; ++j) {
            double x2y2 = x2 + std::pow(y - j, 2);
            for (int k = zm; k < zp; ++k) {
                // fill a sphere of radius [0, vdw] around the atom
                double dist = x2y2 + std::pow(z - k, 2);
                auto& bin = gref.index(i, j, k);
                if (dist <= rvdw2) {
                    if (!gref.is_empty_or_volume(bin)) {continue;}
                    added_volume += !gref.is_volume(bin); // only add to the volume if the bin is not already part of the volume
                    bin = grid::detail::A_AREA;
                }
            }
        }
    }
    grid->add_volume(added_volume);
}

// linear interpolation of the shell width as described in the paper
auto get_shell_width(double Rg) {
    double a = (5. - 3.)/(20.-15.);
    return std::clamp(a*Rg, 3., 5.);
}

#include "data/Body.h"
std::vector<data::record::Water> PepsiHydration::generate_explicit_hydration() {
    double shell_width = get_shell_width(protein->get_Rg());
    double r = 3; // distance from the atom to the hydration shell

    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    std::vector<Water> placed_water;
    placed_water.reserve(grid->a_members.size());
    auto add_loc = [&] (const Vector3<int>& v) {
        placed_water.emplace_back(Water::create_new_water(grid->to_xyz(v)));
    };

    // loop over the location of all member atoms
    double max_r = r+shell_width;
    double max_r2 = max_r*max_r;
    for (const auto& atom : grid->a_members) {
        auto coords_abs = atom.get_atom().get_coordinates();

        // scan for free cells in a box of size [x-r, x+r][y-r, y+r][z-r, z+r]
        auto bin_min = grid->to_bins(coords_abs - max_r);
        auto bin_max = grid->to_bins(coords_abs + max_r);
        bin_min.x() = std::max<int>(bin_min.x()-1, 0); bin_max.x() = std::min<int>(bin_max.x()+1, bins[0]-1);
        bin_min.y() = std::max<int>(bin_min.y()-1, 0); bin_max.y() = std::min<int>(bin_max.y()+1, bins[1]-1);
        bin_min.z() = std::max<int>(bin_min.z()-1, 0); bin_max.z() = std::min<int>(bin_max.z()+1, bins[2]-1);

        for (int i = bin_min.x(); i <= bin_max.x(); ++i) {
            for (int j = bin_min.y(); j <= bin_max.y(); ++j) {
                for (int k = bin_min.z(); k <= bin_max.z(); ++k) {
                    if (!gref.is_empty(i, j, k)) {
                        continue;
                    }

                    double dist = grid->to_xyz(Vector3<int>(i, j, k)).distance2(coords_abs);
                    if (dist < max_r2) {
                        add_loc(Vector3<int>(i, j, k));
                        gref.index(i, j, k) = grid::detail::W_CENTER;
                    }
                }
            }
        }
    }
    return placed_water;
}