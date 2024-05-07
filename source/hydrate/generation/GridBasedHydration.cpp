/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/generation/GridBasedHydration.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <hydrate/ExplicitHydration.h>
#include <hydrate/culling/CullingFactory.h>
#include <data/record/Water.h>
#include <utility/Console.h>
#include <settings/GridSettings.h>

using namespace hydrate;

GridBasedHydration::GridBasedHydration(observer_ptr<grid::Grid> grid) : grid(grid), culling_strategy(hydrate::factory::construct_culling_strategy(grid)) {}
GridBasedHydration::~GridBasedHydration() = default;

std::unique_ptr<Hydration> GridBasedHydration::hydrate() {
    if (grid->w_members.size() != 0) {console::print_warning("Warning in GridBasedHydration::hydrate: Attempting to hydrate a grid which already contains water!");}
    std::vector<data::record::Water> waters = generate_explicit_hydration();

    // assume the protein is a perfect sphere. then we want the number of water molecules to be proportional to the surface area
    double vol = grid->get_volume();                    // volume in cubic Ångström
    double r = std::cbrt(3*vol/(4*constants::pi));      // radius of the protein in Ångström
    double area = 4*constants::pi*std::pow(r, 2.5);     // surface area of the protein in Ångström^2
    double target = settings::grid::water_scaling*area; // the target number of water molecules

    culling_strategy->set_target_count(target);
    waters = culling_strategy->cull(std::move(waters));
    return std::make_unique<ExplicitHydration>(std::move(waters));
}