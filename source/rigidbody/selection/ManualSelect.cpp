/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/ManualSelect.h>
#include <utility/Exceptions.h>

using namespace rigidbody::selection;

ManualSelect::ManualSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {}

ManualSelect::~ManualSelect() = default;

std::pair<unsigned int, unsigned int> ManualSelect::next() {
    throw except::not_implemented("ManualSelect::next: Not implemented.");
}
