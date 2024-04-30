/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/setup/LoadExistingElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody::sequencer;

LoadExistingElement::LoadExistingElement(observer_ptr<Sequencer> owner, observer_ptr<RigidBody> rigidbody) : owner(owner), rigidbody(rigidbody) {
    owner->_set_active_body(rigidbody);
}

void LoadExistingElement::run() {
    owner->_get_rigidbody() = rigidbody;
}