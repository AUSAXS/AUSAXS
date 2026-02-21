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
        assert(index.replica == 0 && index.symmetry == -1 && "SymmetryElement::SymmetryElement: The body name must refer to a base body (symmetry -1, replica 0).");

        // ensure the body's symmetry storage is optimizable
        assert(dynamic_cast<symmetry::OptimizableSymmetryStorage*>(molecule->get_body(ibody).symmetry().get_obj()));
        molecule->get_body(ibody).symmetry().add(symmetry[i]);
        rigidbody->conformation->initial_conformation[ibody].symmetry().add(symmetry[i]);

        // add names for the symmetric bodies
        auto& name_map = setup._get_body_names();
        int isymmetry = molecule->get_body(ibody).size_symmetry()-1;
        assert(0 <= isymmetry && "SymmetryElement::SymmetryElement: Inconsistent data structures.");
        if (int reps = molecule->get_body(ibody).symmetry().get(isymmetry)->repetitions(); reps == 1) { // single replica only: b1s1
            name_map.emplace(names[i] + "s" + std::to_string(isymmetry+1), detail::to_index(ibody, isymmetry, 0));
        } else { // multiple replicas, so include replica index in name: b1s1r1, b1s1r2, ...
            for (int j = 0; j < reps; ++j) {
                name_map.emplace(names[i] + "s" + std::to_string(isymmetry+1) + "r" + std::to_string(j+1), detail::to_index(ibody, isymmetry, j));
            }
        }

        // place the symmetry body at a sane distance from the original
        *molecule->get_body(ibody).symmetry().get(isymmetry)->span_translation().begin() = 2*molecule->get_Rg(false);
        *rigidbody->conformation->initial_conformation[ibody].symmetry().get(isymmetry)->span_translation().begin() = 2*molecule->get_Rg(false);
        rigidbody->conformation->absolute_parameters.parameters[ibody].symmetry_pars.emplace_back(
            molecule->get_body(ibody).symmetry().get(isymmetry)->clone()
        );
    }    
    // Refresh grid to accommodate symmetry bodies
    rigidbody->refresh_grid();
}

SymmetryElement::~SymmetryElement() = default;

void SymmetryElement::run() {}