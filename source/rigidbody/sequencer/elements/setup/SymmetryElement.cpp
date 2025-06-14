// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/SymmetryElement.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/HistogramSettings.h>
#include <data/Body.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

SymmetryElement::SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<symmetry::type>& symmetry) 
    : owner(owner) 
{
    assert(names.size() == symmetry.size() && "SymmetryElement::SymmetryElement: The number of names and symmetries must be equal.");

    owner->_get_molecule()->set_histogram_manager(std::make_unique<hist::PartialSymmetryManagerMT<true>>(owner->_get_molecule()));
    for (unsigned int i = 0; i < names.size(); ++i) {
        if (!owner->setup()->_get_body_names().contains(names[i])) {
            std::cout << "Body names:" << std::endl;
            for (const auto& [name, index] : owner->setup()->_get_body_names()) {
                std::cout << name << " " << index << std::endl;
            }
            throw std::runtime_error("SymmetryElement::SymmetryElement: The body name \"" + names[i] + "\" is not known.");
        }
        int ibody = owner->setup()->_get_body_names().at(names[i]);

        // the body symmetry storage must be replaced with an OptimizableSymmetryStorage object
        if (auto obj = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(owner->_get_molecule()->get_body(ibody).symmetry().get_obj()); !obj) {
            // sanity check: ensure that the body does not already have symmetries before replacing the object
            assert(owner->_get_molecule()->get_body(ibody).size_symmetry() == 0 && "SymmetryElement::SymmetryElement: The body already has symmetries.");
            owner->_get_molecule()->get_body(ibody).symmetry().set_obj(std::make_unique<symmetry::OptimizableSymmetryStorage>());
        }
        owner->_get_molecule()->get_body(ibody).symmetry().add(symmetry[i]);

        // place the symmetry body at a sane distance from the original
        double Rg = owner->_get_molecule()->get_Rg();
        owner->_get_molecule()->get_body(ibody).symmetry().get(0).initial_relation.translation = {2*Rg, 0, 0};

        std::cout << "SymmetryElement::SymmetryElement: Added symmetry to body " << names[i] << std::endl;
        std::cout << "\tIt now has " << owner->_get_molecule()->get_body(ibody).size_symmetry() << " symmetries." << std::endl;
    }
}

SymmetryElement::~SymmetryElement() = default;

void SymmetryElement::run() {}