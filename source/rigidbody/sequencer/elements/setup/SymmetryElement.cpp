// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/Rigidbody.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <rigidbody/sequencer/elements/setup/SymmetryElement.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/HistogramSettings.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

SymmetryElement::SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<symmetry::type>& symmetry) 
    : owner(owner) 
{
    assert(names.size() == symmetry.size() && "SymmetryElement::SymmetryElement: The number of names and symmetries must be equal.");
    auto molecule = owner->_get_molecule();
    auto rigidbody = owner->_get_rigidbody();
    auto& setup = owner->setup();

    molecule->set_histogram_manager(std::make_unique<hist::PartialSymmetryManagerMT<true, false>>(molecule));
    for (unsigned int i = 0; i < names.size(); ++i) {
        auto index = setup._get_body_index(names[i]);
        int ibody = index.body;
        assert(index.replica == 0 && index.symmetry == 0 && "SymmetryElement::SymmetryElement: The body name must refer to a base body (replica 0, symmetry 0).");

        // ensure the body's symmetry storage is optimizable
        assert(dynamic_cast<symmetry::OptimizableSymmetryStorage*>(molecule->get_body(ibody).symmetry().get_obj()));
        molecule->get_body(ibody).symmetry().add(symmetry[i]);
        rigidbody->conformation->initial_conformation[ibody].symmetry().add(symmetry[i]);

        // add names for the symmetry bodies
        auto& name_map = setup._get_body_names();
        int isymmetry = molecule->get_body(ibody).size_symmetry();
        if (int reps = molecule->get_body(ibody).symmetry().get(isymmetry-1).repetitions; reps == 1) { // single replica only: b1s1
            name_map.emplace(names[i] + "s" + std::to_string(isymmetry), detail::to_index(ibody, isymmetry-1, i-1));
        } else { // multiple replicas, so include replica index in name: b1s1r1, b1s1r2, ...
            for (int i = 0; i < reps; ++i) {
                name_map.emplace(names[i] + "s" + std::to_string(isymmetry) + "r" + std::to_string(i+1), detail::to_index(ibody, isymmetry-1, i-1));
            }
        }

        // place the symmetry body at a sane distance from the original
        double Rg = molecule->get_Rg();
        molecule->get_body(ibody).symmetry().get(0).initial_relation.translation = {2*Rg, 0, 0};
        rigidbody->conformation->initial_conformation[ibody].symmetry().get(0).initial_relation.translation = {2*Rg, 0, 0};
        rigidbody->conformation->absolute_parameters.parameters[ibody].symmetry_pars.emplace_back(
            molecule->get_body(ibody).symmetry().get(0)
        );

    }
    
    // Refresh grid to accommodate symmetry bodies
    rigidbody->refresh_grid();
}

SymmetryElement::~SymmetryElement() = default;

void SymmetryElement::run() {}