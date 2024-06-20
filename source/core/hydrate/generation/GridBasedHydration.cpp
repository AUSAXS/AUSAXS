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
#include <data/Molecule.h>
#include <utility/Console.h>
#include <settings/GridSettings.h>

#include <cassert>

using namespace hydrate;

GridBasedHydration::GridBasedHydration(observer_ptr<data::Molecule> protein) : protein(protein), culling_strategy(factory::construct_culling_strategy(protein)) {}
GridBasedHydration::GridBasedHydration(observer_ptr<data::Molecule> protein, std::unique_ptr<CullingStrategy> culling_strategy) : protein(protein), culling_strategy(std::move(culling_strategy)) {}

void GridBasedHydration::initialize() {
    protein->signal_modified_hydration_layer();
    if (auto grid = protein->get_grid(); grid == nullptr) {protein->create_grid();}
    else {grid->clear_waters();}
}

GridBasedHydration::~GridBasedHydration() = default;

void GridBasedHydration::set_culling_strategy(std::unique_ptr<CullingStrategy> culling_strategy) {
    this->culling_strategy = std::move(culling_strategy);
}

std::unique_ptr<Hydration> GridBasedHydration::hydrate() {
    assert(protein != nullptr && "GridBasedHydration::hydrate: protein is nullptr");

    auto grid = protein->get_grid();
    assert(grid != nullptr && "GridBasedHydration::hydrate: grid is nullptr");

    if (grid->w_members.size() != 0) {grid->clear_waters();}
    grid->expand_volume();
    auto waters = generate_explicit_hydration();

    // assume the protein is a perfect sphere. then we want the number of water molecules to be proportional to the surface area
    double vol = grid->get_volume();                    // volume in cubic Ångström
    double r = std::cbrt(3*vol/(4*constants::pi));      // radius of the protein in Ångström
    double area = 4*constants::pi*std::pow(r, 2.5);     // surface area of the protein in Ångström^2
    double target = settings::grid::water_scaling*area; // the target number of water molecules

    culling_strategy->set_target_count(target);
    auto remaining_waters = culling_strategy->cull(waters);
    return std::make_unique<ExplicitHydration>(std::move(remaining_waters));
}