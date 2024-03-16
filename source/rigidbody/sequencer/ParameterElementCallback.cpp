/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/ParameterElementCallback.h>
#include <rigidbody/sequencer/ParameterElement.h>

using namespace rigidbody::sequencer;

ParameterElementCallback::ParameterElementCallback(ParameterElement* caller) : LoopElementCallback(caller->owner), caller(caller) {}

ParameterElementCallback::~ParameterElementCallback() = default;

ParameterElement& ParameterElementCallback::amplitude(double amplitude) {
    return caller->amplitude(amplitude);
}