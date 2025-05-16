/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/generation/PepsiHydration.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/Molecule.h>
#include <settings/GridSettings.h>
#include <settings/MoleculeSettings.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::hydrate;

PepsiHydration::PepsiHydration(observer_ptr<data::Molecule> protein) : GridBasedHydration(protein) {
    initialize();
}

PepsiHydration::PepsiHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy) : GridBasedHydration(protein, std::move(culling_strategy)) {
    initialize();
}

void PepsiHydration::initialize() {
    settings::hydrate::culling_strategy = settings::hydrate::CullingStrategy::NoStrategy;
    GridBasedHydration::initialize();
}

PepsiHydration::~PepsiHydration() = default;

// linear interpolation of the shell width as described in the paper
auto get_shell_width(double Rg) {
    double a = (5. - 3.)/(20.-15.);
    return std::clamp(a*Rg, 3., 5.);
}

std::span<grid::GridMember<data::Water>> PepsiHydration::generate_explicit_hydration(std::span<grid::GridMember<data::AtomFF>> atoms) {
    double shell_width = get_shell_width(protein->get_Rg());
    double r = 3; // distance from the atom to the hydration shell

    assert(protein != nullptr && "PepsiHydration::generate_explicit_hydration: protein is nullptr.");
    auto grid = protein->get_grid();
    assert(grid != nullptr && "PepsiHydration::generate_explicit_hydration: grid is nullptr.");

    grid::detail::GridObj& gref = grid->grid;
    auto bins = grid->get_bins();

    std::vector<data::Water> placed_water;
    placed_water.reserve(atoms.size());
    auto add_loc = [&] (const Vector3<int>& v) {
        placed_water.emplace_back(data::Water(grid->to_xyz(v)));
    };

    // loop over the location of all member atoms
    double max_r = r+shell_width;
    double max_r2 = max_r*max_r;
    for (const auto& atom : atoms) {
        const auto& coords_abs = atom.get_atom().coordinates();

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
                        gref.index(i, j, k) &= grid::detail::W_CENTER;
                    }
                }
            }
        }
    }

    return grid->add(placed_water);
}