// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/generation/GridBasedHydration.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <hydrate/ExplicitHydration.h>
#include <hydrate/culling/CullingFactory.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/Console.h>
#include <settings/GridSettings.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::hydrate;

GridBasedHydration::GridBasedHydration(observer_ptr<data::Molecule> protein) : protein(protein) {}
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

void GridBasedHydration::hydrate() {
    assert(protein != nullptr && "GridBasedHydration::hydrate: protein is nullptr");

    auto grid = protein->get_grid();
    assert(grid != nullptr && "GridBasedHydration::hydrate: grid is nullptr");

    if (!culling_strategy) {culling_strategy = factory::construct_culling_strategy(protein, global());}

    if (grid->w_members.size() != 0) {grid->clear_waters();}
    grid->expand_volume();

    // assume the protein is a perfect sphere. then we want the number of water molecules to be proportional to the surface area
    double vol = grid->get_volume();                    // volume in cubic Ångström
    double r = std::cbrt(3*vol/(4*std::numbers::pi));   // radius of the protein in Ångström
    double area = 4*std::numbers::pi*std::pow(r, 2.5);  // surface area of the protein in Ångström^2
    double target = settings::grid::water_scaling*area; // the target number of water molecules

    auto to_atoms = [] (std::span<grid::GridMember<data::Water>> waters) {
        std::vector<data::Water> remaining_waters(waters.size());
        std::transform(waters.begin(), waters.end(), remaining_waters.begin(), 
            [] (const auto& water) {return water.get_atom();}
        );
        return remaining_waters;
    };

    if (global()) { // global hydration
        auto waters = generate_explicit_hydration(grid->a_members);
        culling_strategy->set_target_count(target);
        culling_strategy->cull(waters);
        protein->get_body(0).set_hydration(std::make_unique<ExplicitHydration>(to_atoms(waters)));
        return;
    }

    int start = 0; // body-specific hydration
    for (int i = 0; i < static_cast<int>(protein->size_body()); i++) {
        assert(grid->body_start.contains(protein->get_body(i).get_uid()) && "GridBasedHydration::hydrate: body_start does not contain body uid");
        int end = start + protein->get_body(i).size_atom();
        assert(end <= static_cast<int>(grid->a_members.size()) && "GridBasedHydration::hydrate: Contained bodies have been modified after being added to the grid.");

        auto atoms = std::span(grid->a_members.begin() + start, grid->a_members.begin() + end);
        auto waters = generate_explicit_hydration(atoms);
        culling_strategy->set_target_count(target);
        culling_strategy->cull(waters);
        protein->get_body(i).set_hydration(std::make_unique<ExplicitHydration>(to_atoms(waters)));
        start = end;
    }
}