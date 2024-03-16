/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/ConstraintIteratorElementCallback.h>
#include <rigidbody/sequencer/ConstraintIteratorElement.h>

using namespace rigidbody::sequencer;

ConstraintIteratorElementCallback::ConstraintIteratorElementCallback(ConstraintIteratorElement* caller) : LoopElementCallback(caller->owner), caller(caller) {}