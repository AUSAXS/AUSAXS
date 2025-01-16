/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/setup/SymmetryElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <data/Body.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

SymmetryElement::SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<symmetry::type>& symmetry) 
    : owner(owner) 
{
    assert(names.size() == symmetry.size() && "SymmetryElement::SymmetryElement: The number of names and symmetries must be equal.");

    for (unsigned int i = 0; i < names.size(); ++i) {
        if (!owner->_get_body_names().contains(names[i])) {
            throw std::runtime_error("SymmetryElement::SymmetryElement: The body name \"" + names[i] + "\" is not known.");
        }
        int ibody = owner->_get_body_names().at(names[i]);

        // the body symmetry storage must be replaced with an OptimizableSymmetryStorage object
        if (auto obj = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(owner->_get_rigidbody()->get_body(ibody).symmetry().get_obj()); !obj) {
            // sanity check: ensure that the body does not already have symmetries before replacing the object
            assert(owner->_get_rigidbody()->get_body(ibody).size_symmetry() == 0 && "SymmetryElement::SymmetryElement: The body already has symmetries.");
            owner->_get_rigidbody()->get_body(ibody).symmetry().set_obj(std::make_unique<symmetry::OptimizableSymmetryStorage>());
        }
        owner->_get_rigidbody()->get_body(ibody).symmetry().add(symmetry[i]);
    }
}

SymmetryElement::~SymmetryElement() = default;

void SymmetryElement::run() {}