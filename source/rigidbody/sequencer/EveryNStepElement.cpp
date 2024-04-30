/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/EveryNStepElement.h>

using namespace rigidbody::sequencer;

EveryNStepElement::EveryNStepElement(observer_ptr<LoopElement> owner, unsigned int n) : LoopElement(owner, 1), n(n), loop_counter(0) {}

EveryNStepElement::~EveryNStepElement() = default;

void EveryNStepElement::run() {
    if (++loop_counter % n == 0) {
        for (auto& e : elements) {
            e->run();
        }
    }
}