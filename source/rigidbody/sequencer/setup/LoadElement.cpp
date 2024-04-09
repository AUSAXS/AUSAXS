/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/setup/LoadElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody::sequencer;

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths) : owner(owner) {
    rigidbody = std::make_unique<RigidBody>(data::Molecule(paths));
}

void LoadElement::run() {
    owner->_get_rigidbody() = rigidbody.get();
}